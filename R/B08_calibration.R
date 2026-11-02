# ---------------------------------------------------------------------------
# B08_calibration.R  -  block B8: history matching against the calibration targets.
#
# Validation targets are not read here, and the runner refuses this block if this file so
# much as names a validation targets file.
#
#   * implausibility I = |model - observed| / sqrt(V_obs + V_md) for each of the 26
#     calibration statistics, with V_md = (0.10 x observed)^2 and V_sim = 0 (the population
#     is integrated on a grid, not sampled); for a lower bound only a shortfall counts;
#   * for each parameter set, every population on the search grid (mu, sigma, sd_tau) whose
#     largest I is at most 3 and second largest at most 2.5, and which meets requirements R1
#     and R3, is kept;
#   * waves continue until at least 500 parameter sets are accepted, or until three
#     successive waves each add fewer than 5% to the accepted set and each lower the best
#     largest implausibility by less than 0.05, with at most 8 waves;
#   * after the first wave, half of each wave is drawn in a small box around EACH of the 30
#     best (or accepted) sets, a quarter in a box three times as wide around each, and a
#     quarter from the whole prior;
#   * C5_p_presented is a sensitivity analysis, recorded but not fitted;
#   * acceptance: for every target, the ensemble's median prediction has I at most 2.  If
#     not, the analysis stops here and that is reported.
#
# The ensemble has one member per accepted parameter set, with a population drawn at random
# from the populations kept for it.  B08_pairs.csv also keeps, flagged, the populations that
# pass with targets left out or with alternative readings of C1, C2 and C5.
#
# Outputs: output/B08_targets.csv, output/B08_all.csv, output/B08_waves.csv,
#          output/B08_pairs.csv, output/B08_fine_grid.csv, output/B08_c6_mix.csv,
#          output/B08_ensemble.csv, output/B08_acceptance.csv, output/B08_intervals.csv,
#          output/B08_calibration.json
# ---------------------------------------------------------------------------

source(file.path(QSP2_DIR, "R", "qsp2_model.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_lhs.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_pop.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_cal.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_parallel.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "run_tests.R"), local = TRUE)
build_qsp2_dll(); load_qsp2_dll()

cat("B8  calibration by history matching\n\n")
tr <- qsp2_run_tests()
check("the model still passes its regression tests", isTRUE(tr$ok),
      sprintf("%d tests, %d failed, %d errors", tr$n, tr$failed, tr$errors))

targets <- read_targets("calibration")
check("the calibration targets are available", nrow(targets) == 6,
      paste(targets$id, collapse = ", "))

# C1's variances come from a bootstrap of the case file; only cases with both dates complete
# to the day are used
c1 <- qsp2_c1_bootstrap(QSP2_C1_FILE())
check("the case file gives the C1 quartiles",
      c1$complete && c1$n_mono == 988 && c1$n_combo == 262 &&
        max(abs(c1$q_combo - c(14.25, 20, 39.75))) < 1e-9 && max(abs(c1$q_mono - c(20, 30, 89))) < 1e-9,
      sprintf("combination %s (n = %d of %d with any date), monotherapy %s (n = %d of %d); bootstrap SD of the ratio %.4f",
              paste(c1$q_combo, collapse = "/"), c1$n_combo, c1$n_combo_any_date,
              paste(c1$q_mono, collapse = "/"), c1$n_mono, c1$n_mono_any_date, c1$sd_ratio))
# the horizon C1 is computed over must be the one the model simulates it over, so the two
# cannot drift apart
check("C1 is computed over the horizon cohort S0 is simulated over",
      identical(as.numeric(c1$horizon), as.numeric(qsp2_cal_cohorts()$S0$H)),
      sprintf("%s days; %d combination and %d monotherapy onsets fall after it and are left out",
              c1$horizon, c1$beyond_combo, c1$beyond_mono))
TGc1 <- qsp2_cal_targets(c1)
check("the C1 targets are the case file's own values within the horizon",
      max(abs(TGc1$obs[match(c("C1_q25_combo", "C1_q50_combo", "C1_q75_combo"), TGc1$id)] - c1$q_combo)) < 1e-9 &&
        max(abs(TGc1$obs[match(c("C1_q25_mono", "C1_q50_mono", "C1_q75_mono"), TGc1$id)] - c1$q_mono)) < 1e-9 &&
        abs(TGc1$obs[TGc1$id == "C1_ratio"] - c1$q_mono[2] / c1$q_combo[2]) < 1e-9,
      sprintf("ratio of medians %.4f", TGc1$obs[TGc1$id == "C1_ratio"]))
TG <- qsp2_cal_targets(c1)
TGS <- qsp2_cal_targets_sens()
fwrite(rbind(cbind(TG, use = "calibration"), cbind(TGS, use = "sensitivity only")),
       file.path(OUT_DIR, "B08_targets.csv"))
check("every target has an observed value and a positive total variance",
      all(is.finite(TG$obs)) && all(TG$v_obs + TG$v_md > 0),
      sprintf("%d statistics: %s", nrow(TG), paste(TG$id, collapse = ", ")))
check("the C5 grade target is a sensitivity analysis, not a calibration target",
      !("C5_p_presented" %in% TG$id) && identical(TGS$id, "C5_p_presented") && nrow(TG) == 26,
      sprintf("%d calibration targets; recorded only: %s", nrow(TG), paste(TGS$id, collapse = ", ")))

feas <- fread(file.path(OUT_DIR, "B06_feasible.csv"))
check("the region B6 found is available to start from", nrow(feas) >= 10, sprintf("%d sets", nrow(feas)))

CTX <- list(mode = "B8", cohorts = qsp2_cal_cohorts(), sgrid = QSP2_SGRID_CAL,
            tri = QSP2_TRI_GRID, tg = TG, tg_sens = TGS, ipi = FALSE)
NPW  <- as.integer(Sys.getenv("QSP2_B8_N", "600"))
GOAL <- 500L
MAXW <- 8L
DEL_E <- qsp2_default_pars()[["del_E"]]
valid <- function(d) d[rho < DEL_E + k_tr]

run_wave <- function(draws, w, cl) {
  jobs <- lapply(seq_len(nrow(draws)), function(i) list(id = i, draw = as.list(draws[i])))
  t0 <- Sys.time()
  res <- qsp2_papply(cl, jobs, "qsp2_cal_one", CTX)
  st <- vapply(res, function(x) x$status, "")
  ok <- which(st == "ok")
  best <- rbindlist(lapply(ok, function(i) as.data.table(as.list(res[[i]]$best))), fill = TRUE)
  best <- cbind(set = paste0(w, "_", ok), wave = w, draws[ok], best)
  kept <- rbindlist(lapply(ok, function(i) {
    k <- res[[i]]$keep
    if (is.null(k)) return(NULL)
    cbind(set = paste0(w, "_", i), wave = w, draws[rep(i, nrow(k))], as.data.table(k))
  }), fill = TRUE)
  acc_rows <- if (nrow(kept)) kept[pass_all == 1] else kept
  cat(sprintf("wave %d: %d sets, %d ran, %d errors; %d sets accepted with %d populations; best largest I %.2f (%.0f min)\n",
              w, nrow(draws), length(ok), sum(grepl("^error", st)),
              if (nrow(acc_rows)) uniqueN(acc_rows$set) else 0L, nrow(acc_rows),
              suppressWarnings(min(best$Imax, na.rm = TRUE)),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  errs <- unique(st[grepl("^error", st)])
  if (length(errs)) cat("  errors:", paste(head(errs, 3), collapse = " | "), "\n")
  list(best = best, kept = kept, acc = acc_rows, n_err = sum(grepl("^error", st)))
}

set.seed(20260929)
cl <- qsp2_cluster()
all_best <- NULL; pairs <- NULL; wlog <- NULL
PRI <- names(QSP2_PRIOR)
small <- 0L; w <- 0L
repeat {
  w <- w + 1L
  n_acc <- if (is.null(pairs) || !nrow(pairs)) 0L else uniqueN(pairs[pass_all == 1]$set)
  draws <- if (w == 1L) {
    # four fifths from the region B6 found, jittered; one fifth from the whole prior
    nh <- as.integer(round(0.8 * NPW))
    d1 <- feas[sample.int(nrow(feas), 4L * nh, replace = TRUE), ..PRI]
    for (nm in PRI) {
      spec <- QSP2_PRIOR[[nm]]; lo <- as.numeric(spec[1]); hi <- as.numeric(spec[2])
      v <- d1[[nm]]
      v <- if (spec[3] == "log") exp(log(v) + rnorm(length(v), 0, 0.25)) else
        v + rnorm(length(v), 0, 0.05 * (hi - lo))
      if (spec[3] == "int") v <- round(v)
      set(d1, j = nm, value = pmin(pmax(v, lo), hi))
    }
    d1 <- valid(d1)[seq_len(nh)]
    d2 <- valid(qsp2_lhs(4L * (NPW - nh), QSP2_PRIOR))[seq_len(NPW - nh)]
    rbind(d1, d2)
  } else {
    # around the accepted sets when there are enough of them, otherwise around the 30 sets
    # with the lowest largest implausibility: a box around EACH point (half the wave at 10% of
    # each prior range either side, a quarter at 30%) and a quarter from the whole prior.
    around <- if (n_acc >= 10) unique(pairs[pass_all == 1, c("set", PRI), with = FALSE])[, ..PRI] else
      all_best[order(Imax)][1:min(30, nrow(all_best)), ..PRI]
    n1 <- NPW %/% 2L; n2 <- NPW %/% 4L; n3 <- NPW - n1 - n2
    rbind(valid(qsp2_lhs_near(4L * n1, QSP2_PRIOR, around, 0.1))[seq_len(n1)],
          valid(qsp2_lhs_near(4L * n2, QSP2_PRIOR, around, 0.3))[seq_len(n2)],
          valid(qsp2_lhs(4L * n3, QSP2_PRIOR))[seq_len(n3)])
  }
  draws <- draws[complete.cases(draws)]
  best_before <- if (is.null(all_best)) Inf else suppressWarnings(min(all_best$Imax, na.rm = TRUE))
  r <- run_wave(draws, w, cl)
  all_best <- rbind(all_best, r$best, fill = TRUE)
  pairs <- rbind(pairs, r$kept, fill = TRUE)
  n_new <- if (nrow(r$acc)) uniqueN(r$acc$set) else 0L
  growth <- n_new / max(n_acc, 1)
  best_now <- suppressWarnings(min(all_best$Imax, na.rm = TRUE))
  gain <- if (is.finite(best_before)) best_before - best_now else Inf
  # a wave counts as small only if it adds fewer than 5% AND lowers the best largest
  # implausibility by less than 0.05, so a search that is still improving goes on
  small <- if (growth < 0.05 && gain < 0.05) small + 1L else 0L
  wlog <- rbind(wlog, data.table(wave = w, drawn = nrow(draws), ran = nrow(r$best), errors = r$n_err,
                                 accepted_sets = n_new, accepted_populations = nrow(r$acc),
                                 total_accepted = n_acc + n_new, growth = growth,
                                 best_Imax = suppressWarnings(min(r$best$Imax, na.rm = TRUE)),
                                 best_so_far = best_now, gain = gain))
  if (n_acc + n_new >= GOAL) { stop_reason <- "at least 500 sets accepted"; break }
  if (small >= 3L) { stop_reason <- "three successive waves each added fewer than 5% and improved the best largest implausibility by less than 0.05"; break }
  if (w >= MAXW) { stop_reason <- sprintf("the limit of %d waves", MAXW); break }
}
# the populations accepted under the acceptance rule (the pool also holds those kept only for
# the alternative readings); an empty pool has no columns at all
accepted <- function() if (is.null(pairs) || !nrow(pairs) || !"pass_all" %in% names(pairs))
  data.table() else pairs[pass_all == 1]

# Acceptance is decided on a fine S grid.  On two test sets the 25-node calibration grid and a
# 41-node grid differed by up to 0.23 in the second-largest
# implausibility of the best population, and a 41-node and an 81-node grid by up to 0.08, so a set
# near the cut-off could be accepted or rejected by the grid alone.  The waves search on the
# 25-node grid; every set they accept, and every set whose best population comes within a margin
# wider than that error (largest I at most 3.5, second largest at most 3.0), is evaluated again on
# 81 nodes, and only those results decide acceptance and make the ensemble.
SG_FINE <- 10^seq(-4, 4, length.out = 81)
acc_coarse <- if (nrow(accepted())) unique(accepted()$set) else character(0)
cand <- all_best[set %in% union(acc_coarse, all_best[is.finite(Imax) & Imax <= 3.5 & I2nd <= 3.0, set])]
fine_log <- data.table()
if (nrow(cand)) {
  CTXF <- CTX; CTXF$sgrid <- SG_FINE
  t0f <- Sys.time()
  rf <- qsp2_papply(cl, lapply(seq_len(nrow(cand)), function(i)
    list(id = i, draw = as.list(cand[i, ..PRI]))), "qsp2_cal_one", CTXF)
  okf <- vapply(rf, function(x) identical(x$status, "ok"), TRUE)
  keptF <- rbindlist(lapply(which(okf), function(i) {
    k <- rf[[i]]$keep
    if (is.null(k)) return(NULL)
    cbind(set = cand$set[i], wave = cand$wave[i], cand[rep(i, nrow(k)), ..PRI], as.data.table(k))
  }), fill = TRUE)
  acc_fine <- if (nrow(keptF) && "pass_all" %in% names(keptF)) unique(keptF[pass_all == 1]$set) else character(0)
  bestF <- function(nm) vapply(seq_along(rf), function(i) if (okf[i]) unname(rf[[i]]$best[[nm]]) else NA_real_, 0)
  fine_log <- data.table(set = cand$set, wave = cand$wave, Imax_25 = cand$Imax, I2nd_25 = cand$I2nd,
                         Imax_81 = bestF("Imax"), I2nd_81 = bestF("I2nd"),
                         accepted_25 = cand$set %in% acc_coarse, accepted_81 = cand$set %in% acc_fine)
  fwrite(fine_log, file.path(OUT_DIR, "B08_fine_grid.csv"))
  # the candidates' populations from the coarse grid are replaced by those from the fine grid
  pairs <- if (is.null(pairs) || !nrow(pairs)) keptF else rbind(pairs[!(set %in% cand$set)], keptF, fill = TRUE)
  cat(sprintf("\nfine grid (81 S nodes): %d candidate sets, %d accepted on 25 nodes, %d on 81 (%.0f min)\n",
              nrow(cand), length(acc_coarse), length(acc_fine),
              as.numeric(difftime(Sys.time(), t0f, units = "mins"))))
}

# C6 (Pereyra Pietri) does not report its regimen mix and 10% combination is assumed; how far
# the statistic moves with the mix, on up to 30 accepted members
MIX <- c(0, 0.1, 0.3)
tropmix <- NULL
if (nrow(accepted())) {
  one <- accepted()[, .SD[1], by = set]
  one <- one[sample.int(nrow(one), min(30, nrow(one)))]
  tm <- qsp2_papply(cl, lapply(seq_len(nrow(one)), function(i)
    list(draw = as.list(one[i, ..PRI]), tri = as.data.frame(one[i, .(mu, sigma, sd_tau)]))),
    "qsp2_cal_tropmix", list(sgrid = SG_FINE, mix = MIX))
  tropmix <- do.call(rbind, tm)
  colnames(tropmix) <- paste0("combo_", MIX)
}
qsp2_cluster_stop(cl)
cat(sprintf("\nwaves stopped: %s\n", stop_reason))
if (!is.null(tropmix)) {
  rel <- (tropmix[, 3] - tropmix[, 1]) / tropmix[, 2]
  cat(sprintf("C6 troponin T at diagnosis: median %.0f ng/L at 0%% combination, %.0f at 10%%, %.0f at 30%%; the change from 0%% to 30%% is %.1f%% of the 10%% value (median over %d members)\n",
              median(tropmix[, 1], na.rm = TRUE), median(tropmix[, 2], na.rm = TRUE),
              median(tropmix[, 3], na.rm = TRUE), 100 * median(rel, na.rm = TRUE), nrow(tropmix)))
  fwrite(as.data.table(tropmix), file.path(OUT_DIR, "B08_c6_mix.csv"))
}
fwrite(all_best, file.path(OUT_DIR, "B08_all.csv"))
fwrite(wlog, file.path(OUT_DIR, "B08_waves.csv"))
if (!is.null(pairs) && nrow(pairs)) fwrite(pairs, file.path(OUT_DIR, "B08_pairs.csv"))

n_sets <- if (nrow(accepted())) uniqueN(accepted()$set) else 0L
check("the waves ran to the stopping rule", TRUE,
      sprintf("%d waves, %d sets evaluated; %s", nrow(wlog), nrow(all_best), stop_reason))
check("almost every set ran", sum(wlog$errors) <= 0.01 * sum(wlog$ran),
      sprintf("%d errors in %d runs", sum(wlog$errors), sum(wlog$ran)))
check("acceptance was decided on the fine S grid",
      !length(acc_coarse) || (nrow(fine_log) > 0 && all(acc_coarse %in% fine_log$set)),
      if (nrow(fine_log)) sprintf("%d candidate sets re-evaluated on 81 nodes: %d accepted on 25 nodes, %d on 81",
                                  nrow(fine_log), sum(fine_log$accepted_25), sum(fine_log$accepted_81))
      else "no set came near acceptance, so none needed re-evaluating")
check("a non-implausible ensemble was found", n_sets > 0,
      sprintf("%d sets with %d populations", n_sets, nrow(accepted())))
# the gate: at least 500 non-implausible parameter sets
check("the ensemble has at least 500 non-implausible sets", n_sets >= 500,
      sprintf("%d sets (at least 500 required)", n_sets))

if (n_sets > 0) {
  # one member per accepted set, its population drawn from those kept for it
  ens <- accepted()[, .SD[sample.int(.N, 1)], by = set]
  ens[, member := seq_len(.N)]
  fwrite(ens, file.path(OUT_DIR, "B08_ensemble.csv"))
  check("every member meets requirements R1 and R3",
        all(qsp2_req_R1(ens)) && all(qsp2_req_R3(ens)), sprintf("%d members", nrow(ens)))

  # acceptance: the ensemble's median prediction, target by target
  med <- vapply(TG$id, function(id) median(ens[[id]], na.rm = TRUE), 0)
  Imed <- qsp2_cal_I(matrix(med, 1, dimnames = list(NULL, TG$id)), TG)[1, ]
  acc <- data.table(target = TG$id, observed = TG$obs, kind = TG$kind, ensemble_median = med,
                    q025 = vapply(TG$id, function(id) quantile(ens[[id]], 0.025, na.rm = TRUE), 0),
                    q975 = vapply(TG$id, function(id) quantile(ens[[id]], 0.975, na.rm = TRUE), 0),
                    I_of_median = Imed)
  fwrite(acc, file.path(OUT_DIR, "B08_acceptance.csv"))
  cat("\nthe ensemble against each target (median and 95% range over the members):\n")
  for (i in seq_len(nrow(acc))) cat(sprintf("  %-16s observed %9.4g   model %9.4g  (%9.4g to %9.4g)   I %5.2f\n",
                                            acc$target[i], acc$observed[i], acc$ensemble_median[i],
                                            acc$q025[i], acc$q975[i], acc$I_of_median[i]))
  check("every target's ensemble median prediction has an implausibility of at most 2",
        all(is.finite(Imed) & Imed <= 2),
        sprintf("largest %.2f, for %s", max(Imed), TG$id[which.max(Imed)]))

  # every fitted quantity with its interval
  PAR_ALL <- c(PRI, "mu", "sigma", "sd_tau")
  intervals <- rbindlist(lapply(PAR_ALL, function(pn) data.table(parameter = pn,
    median = median(ens[[pn]]), lo = quantile(ens[[pn]], 0.025, names = FALSE),
    hi = quantile(ens[[pn]], 0.975, names = FALSE))))
  fwrite(intervals, file.path(OUT_DIR, "B08_intervals.csv"))
  check("every fitted parameter is reported with an interval", all(is.finite(intervals$lo) & is.finite(intervals$hi)),
        sprintf("%d parameters", nrow(intervals)))
}

save_json(list(block = "B08",
               evaluated = nrow(all_best), accepted_sets = n_sets,
               accepted_populations = nrow(accepted()),
               pool_populations = if (is.null(pairs)) 0L else nrow(pairs), stop_reason = stop_reason,
               waves = as.list(as.data.frame(wlog)), targets = TG$id,
               c1_bootstrap_sd = list(combo = c1$sd_combo, mono = c1$sd_mono, ratio = c1$sd_ratio)),
          "B08_calibration.json")
