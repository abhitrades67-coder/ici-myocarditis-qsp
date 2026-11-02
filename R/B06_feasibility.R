# ---------------------------------------------------------------------------
# B06_feasibility.R  -  block B6: can this structure meet requirements R1 to R5?
#
#   R1  among T2 cases, the first T1 is at least 7 days after the first dose in at least 95%
#       of cases (the cohorts' times to detection are fitted directly in block B8);
#   R2  with one parameter set, the population reproduces the incidence values of C2 to C5
#       within the implausibility cut-offs (largest at most 3, second largest at most 2.5).
#       The population of each set is the one on the search grid with the lowest largest
#       implausibility over those values;
#   R3  on that population, with diagnosis at clinical presentation, fewer than half of T2
#       cases reach T3.  The nadir-LVEF bands of severe cases are reported below as an output
#       of the model;
#   R4  somewhere in the feasible region, T3 cases present earlier than T2 cases without T3;
#   R5  with ipilimumab alone, T2 incidence is above zero somewhere in the prior.
# R1 and R3 are read on monotherapy, because the cohorts that set R1's range are 87% to 98%
# monotherapy; the combination values are kept in the output.
#
# The first Latin hypercube has 2,000 sets over the whole prior.  If too few sets meet R1 to
# R3, further waves are drawn around the sets that came closest.  If no region meets R1 to
# R5, the analysis stops here and that is the result.
#
# Outputs: output/B06_lhs.csv, output/B06_feasible.csv, output/B06_feasibility.json
# ---------------------------------------------------------------------------

source(file.path(QSP2_DIR, "R", "qsp2_model.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_lhs.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_pop.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_cal.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_parallel.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "run_tests.R"), local = TRUE)
build_qsp2_dll(); load_qsp2_dll()

cat("B6  structural feasibility: requirements R1 to R5\n\n")
tr <- qsp2_run_tests()
check("the model still passes its regression tests", isTRUE(tr$ok),
      sprintf("%d tests, %d failed, %d errors", tr$n, tr$failed, tr$errors))

N0  <- as.integer(Sys.getenv("QSP2_B6_N0", "2000"))     # the first hypercube
NW  <- as.integer(Sys.getenv("QSP2_B6_N", "1000"))      # each later wave
MAXW <- as.integer(Sys.getenv("QSP2_B6_WAVES", "4"))
ENOUGH <- 30L                                           # sets meeting R1 to R3 to stop early

c1 <- qsp2_c1_bootstrap(QSP2_C1_FILE(), B = 2000)
TG <- qsp2_cal_targets(c1)
R2_IDS <- TG$id[TG$r2]
cat("R2 uses:", paste(R2_IDS, collapse = ", "), "\n\n")
# R2 needs the incidence values of C2 to C5, and R1, R3 and R4 the clinical course (S0), so the
# troponin T cohorts of C6 are not run here
CTX <- list(mode = "B6", cohorts = qsp2_cal_cohorts()[c("S0", "C3", "C4", "C5")], sgrid = QSP2_SGRID_CAL,
            tri = QSP2_TRI_GRID, tg = TG, ipi = FALSE)

#' How far a set is from meeting R1 to R3, on a continuous scale; zero means it meets them.
distance <- function(r) with(r, {
  d2 <- pmax(0, I2max - 3) + pmax(0, I2second - 2.5)
  d1 <- 10 * pmax(0, 0.95 - ifelse(is.finite(R1_p7_mono), R1_p7_mono, 0))
  sh <- ifelse(is.finite(R3_t3_share_mono), R3_t3_share_mono, 1)
  d3 <- 5 * pmax(0, sh - 0.5)
  d1 + d2 + d3
})

run_wave <- function(draws, w, cl) {
  jobs <- lapply(seq_len(nrow(draws)), function(i) list(id = i, draw = as.list(draws[i])))
  t0 <- Sys.time()
  res <- qsp2_papply(cl, jobs, "qsp2_cal_one", CTX)
  st <- vapply(res, function(x) x$status, "")
  ok <- which(st == "ok")
  tab <- rbindlist(lapply(ok, function(i) as.data.table(as.list(res[[i]]$best))), fill = TRUE)
  tab <- cbind(draws[ok], tab, secs = vapply(res[ok], function(x) x$secs, 0),
               failed_nodes = vapply(res[ok], function(x) x$failed_nodes, 0L))
  tab[, wave := w]
  tab[, R1 := qsp2_req_R1(tab)]
  tab[, R2 := is.finite(I2max) & I2max <= 3 & I2second <= 2.5]
  tab[, R3 := qsp2_req_R3(tab)]
  tab[, R4 := qsp2_req_R4(tab)]
  tab[, dist := distance(tab)]
  cat(sprintf("wave %d: %d sets drawn, %d ran, %d broke the rho constraint, %d errors (%.0f min)\n",
              w, nrow(draws), length(ok), sum(st == "constraint"), sum(grepl("^error", st)),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  errs <- unique(st[grepl("^error", st)])
  if (length(errs)) cat("  errors:", paste(head(errs, 3), collapse = " | "), "\n")
  cat(sprintf("  R1 %d, R2 %d, R3 %d, R1 to R3 together %d\n", sum(tab$R1), sum(tab$R2),
              sum(tab$R3), sum(tab$R1 & tab$R2 & tab$R3)))
  list(tab = tab, n_err = sum(grepl("^error", st)), n_con = sum(st == "constraint"))
}

# rho must stay below del_E + k_tr (Section 5 of the specification).  About two draws in three
# from the prior break it, so each wave draws four times what it needs and keeps the first
# sets that respect it; the hypercube's rows are in random order, so that is a random subset.
DEL_E <- qsp2_default_pars()[["del_E"]]
valid <- function(n, sampler) {
  d <- sampler(4L * n)
  d <- d[rho < DEL_E + k_tr]
  d[seq_len(min(n, nrow(d)))]
}

set.seed(20260928)
cl <- qsp2_cluster()
all_res <- NULL; wlog <- NULL
for (w in seq_len(MAXW)) {
  draws <- if (w == 1L) valid(N0, function(n) qsp2_lhs(n, QSP2_PRIOR)) else {
    keep <- all_res[order(dist)][1:min(30, nrow(all_res)), names(QSP2_PRIOR), with = FALSE]
    valid(NW, function(n) qsp2_lhs_in_region(n, QSP2_PRIOR, keep, 0.15))
  }
  r <- run_wave(draws, w, cl)
  all_res <- rbind(all_res, r$tab, fill = TRUE)
  wlog <- rbind(wlog, data.table(wave = w, drawn = nrow(draws), ran = nrow(r$tab),
                                 constraint = r$n_con, errors = r$n_err,
                                 R1 = sum(r$tab$R1), R2 = sum(r$tab$R2), R3 = sum(r$tab$R3),
                                 R1to3 = sum(r$tab$R1 & r$tab$R2 & r$tab$R3),
                                 best_dist = min(r$tab$dist)))
  if (sum(all_res$R1 & all_res$R2 & all_res$R3) >= ENOUGH) break
}

# R5: ipilimumab alone, on the sets that meet R1 to R3 (up to 50) and on 50 sets of the
# first hypercube, which covers the whole prior
feas <- all_res[R1 & R2 & R3]
r5_sets <- rbind(feas[seq_len(min(50, nrow(feas)))],
                 all_res[wave == 1][sample.int(sum(all_res$wave == 1), min(50, sum(all_res$wave == 1)))])
CTX5 <- CTX; CTX5$ipi <- TRUE
r5 <- qsp2_papply(cl, lapply(seq_len(nrow(r5_sets)), function(i)
  list(id = i, draw = as.list(r5_sets[i, names(QSP2_PRIOR), with = FALSE]))), "qsp2_cal_one", CTX5)
qsp2_cluster_stop(cl)
r5_t2 <- vapply(r5, function(x) if (x$status == "ok") unname(x$best[["R5_t2_ipi"]]) else NA_real_, 0)
r5_feas <- r5_t2[seq_len(min(50, nrow(feas)))]

fwrite(all_res, file.path(OUT_DIR, "B06_lhs.csv"))
fwrite(feas, file.path(OUT_DIR, "B06_feasible.csv"))
fwrite(wlog, file.path(OUT_DIR, "B06_waves.csv"))

rates <- c(R1 = mean(all_res$R1), R2 = mean(all_res$R2), R3 = mean(all_res$R3), R4 = mean(all_res$R4))
cat("\nover every wave, share of sets meeting each requirement:\n")
for (nm in names(rates)) cat(sprintf("  %-3s %6.2f%%\n", nm, 100 * rates[[nm]]))
cat(sprintf("  R5  ipilimumab alone gives T2 in %d of %d sets tried (%d of %d that meet R1 to R3)\n",
            sum(r5_t2 > 0, na.rm = TRUE), sum(is.finite(r5_t2)),
            sum(r5_feas > 0, na.rm = TRUE), sum(is.finite(r5_feas))))

check("the first Latin hypercube had at least 2,000 sets over the whole prior",
      wlog$drawn[1] >= 2000,
      sprintf("%d sets respecting the rho constraint in the first wave, %d waves in all, %d sets ran",
              wlog$drawn[1], nrow(wlog), nrow(all_res)))
check("almost every set that respects the rho constraint ran",
      sum(wlog$errors) <= 0.01 * sum(wlog$ran), sprintf("%d errors in %d runs", sum(wlog$errors), sum(wlog$ran)))
for (rq in c("R1", "R2", "R3"))
  check(sprintf("requirement %s is met somewhere in the prior", rq), any(all_res[[rq]]),
        sprintf("%.2f%% of sets", 100 * rates[[rq]]))
check("a region of the prior meets R1, R2 and R3 at the same time", nrow(feas) > 0,
      sprintf("%d sets", nrow(feas)))
check("R4: somewhere in that region, T3 cases present earlier than other T2 cases",
      any(feas$R4), sprintf("%d of %d sets", sum(feas$R4), nrow(feas)))
check("R5: ipilimumab alone gives myocarditis somewhere in the prior",
      any(r5_t2 > 0, na.rm = TRUE),
      sprintf("in %d of %d sets tried, %d of them in the region meeting R1 to R3",
              sum(r5_t2 > 0, na.rm = TRUE), sum(is.finite(r5_t2)), sum(r5_feas > 0, na.rm = TRUE)))
check("the region is large enough to start the calibration from", nrow(feas) >= 10,
      sprintf("%d sets", nrow(feas)))
# the bands are reported, not tested
bandcols <- c("R3_b40_mono", "R3_b30_mono", "R3_b0_mono")
check("the nadir-LVEF bands are reported as an output",
      all(bandcols %in% names(all_res)) && (!nrow(feas) || any(is.finite(unlist(feas[, ..bandcols])))),
      if (nrow(feas)) sprintf("among severe cases in the feasible region: %.0f%% at 40 to 49%%, %.0f%% at 30 to 39%%, %.0f%% below 30%% (medians); no calibration source constrains them",
                              100 * median(feas$R3_b40_mono, na.rm = TRUE),
                              100 * median(feas$R3_b30_mono, na.rm = TRUE),
                              100 * median(feas$R3_b0_mono, na.rm = TRUE))
      else "no feasible sets to report bands for")

if (nrow(feas)) {
  cat("\nin the region meeting R1 to R3 (median and range over the sets):\n")
  show <- c("R1_med_mono", "R1_p7_mono", "R1_med_combo", "R3_t3_share_mono", "R3_b40_mono",
            "R3_b30_mono", "R3_b0_mono", "I2max", "C2_p_combo", "C2_p_mono", "C3_p_damage",
            "C4_p_myo", "C5_p_pos", "C5_ppv", "phi_inf", "mu", "sigma", "sd_tau")
  # a set with no severe case has no LVEF bands, so the summaries skip missing values
  for (nm in show) cat(sprintf("  %-18s %10.4g   (%.4g to %.4g)   %d sets\n", nm,
                               median(feas[[nm]], na.rm = TRUE), min(feas[[nm]], na.rm = TRUE),
                               max(feas[[nm]], na.rm = TRUE), sum(is.finite(feas[[nm]]))))
}

save_json(list(block = "B06",
               first_hypercube = N0, waves = as.list(as.data.frame(wlog)), scored = nrow(all_res),
               rates = as.list(rates), feasible = nrow(feas), feasible_R4 = sum(feas$R4),
               R5_sets_tried = sum(is.finite(r5_t2)), R5_sets_with_T2 = sum(r5_t2 > 0, na.rm = TRUE),
               R2_targets = R2_IDS), "B06_feasibility.json")
