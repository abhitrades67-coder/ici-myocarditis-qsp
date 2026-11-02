# ---------------------------------------------------------------------------
# B07_sensitivity.R  -  block B7: which parameters matter, and which can be estimated.
#
# Two questions, on the sets block B6 scored:
#   1. Sensitivity: which parameters move the calibration statistics?  Partial rank
#      correlation (PRCC) over every set B6 ran, each at the population B6 chose for it.
#   2. Identifiability: how much of each prior survives in the region that meets R1 to R3?
#
# None of the fitted parameters has a literature value; each has a prior range because no
# source gives it.  So none is fixed: every one is carried
# through the calibration with its prior range, and those the data do not pin down are
# reported as such, so no interval in the paper claims more precision than the data give.
#
# Outputs: output/B07_sensitivity.csv, output/B07_profiles.csv, output/B07_sensitivity.json
# ---------------------------------------------------------------------------

source(file.path(QSP2_DIR, "R", "qsp2_model.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_lhs.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_pop.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_cal.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "run_tests.R"), local = TRUE)
build_qsp2_dll(); load_qsp2_dll()

cat("B7  sensitivity and identifiability\n\n")
tr <- qsp2_run_tests()
check("the model still passes its regression tests", isTRUE(tr$ok),
      sprintf("%d tests, %d failed, %d errors", tr$n, tr$failed, tr$errors))

all6 <- fread(file.path(OUT_DIR, "B06_lhs.csv"))
feas <- fread(file.path(OUT_DIR, "B06_feasible.csv"))
check("every set B6 scored is available", nrow(all6) >= 500,
      sprintf("%d sets, of which %d meet R1 to R3", nrow(all6), nrow(feas)))
check("the region meeting R1 to R3 is available", nrow(feas) >= 10, sprintf("%d sets", nrow(feas)))

PARS <- c(names(QSP2_PRIOR), "mu", "sigma", "sd_tau")
OUTS <- c("C1_q50_combo", "C1_q50_mono", "C1_ratio", "C2_p_combo", "C2_p_mono", "C3_p_damage",
          "C3_t_median", "C4_p_myo", "C4_t_mean", "C4_p_screen", "C4_ppv", "C4_p_baseline",
          "C5_p_pos", "C5_q50", "C5_ppv", "C5_t_myo", "C5_p_presented", "R1_med_mono",
          "R3_t3_share_mono")
check("B6 recorded every parameter and every statistic this block needs",
      all(c(PARS, OUTS) %in% names(all6)),
      paste(setdiff(c(PARS, OUTS), names(all6)), collapse = ", "))

prcc <- function(y, x, others) {
  ok <- is.finite(y) & is.finite(x) & Reduce(`&`, lapply(others, is.finite))
  if (sum(ok) < 50) return(NA_real_)
  ry <- rank(y[ok]); rx <- rank(x[ok])
  Z  <- do.call(cbind, lapply(others, function(v) rank(v[ok])))
  suppressWarnings(cor(residuals(lm(ry ~ Z)), residuals(lm(rx ~ Z))))
}
sens <- rbindlist(lapply(OUTS, function(o) rbindlist(lapply(PARS, function(pn)
  data.table(output = o, parameter = pn, n = sum(is.finite(all6[[o]])),
             prcc = prcc(all6[[o]], all6[[pn]], lapply(setdiff(PARS, pn), function(q) all6[[q]])))))))
sens[, abs_prcc := abs(prcc)]
fwrite(sens, file.path(OUT_DIR, "B07_sensitivity.csv"))

infl <- sens[, .(max_abs = suppressWarnings(max(abs_prcc, na.rm = TRUE)),
                 output = output[which.max(abs_prcc)]), by = parameter][order(-max_abs)]
INFLUENTIAL <- infl[is.finite(max_abs) & max_abs >= 0.2]$parameter
NEGLIGIBLE  <- setdiff(PARS, INFLUENTIAL)
cat("\nparameters that move at least one calibration statistic (|PRCC| at least 0.2):\n")
print(infl[parameter %in% INFLUENTIAL][, .(parameter, max_abs = round(max_abs, 3), output)])
cat(sprintf("\n%d move none of them: %s\n", length(NEGLIGIBLE), paste(NEGLIGIBLE, collapse = ", ")))

check("the sensitivity analysis ran for every parameter and statistic",
      nrow(sens) == length(OUTS) * length(PARS) && all(is.finite(sens$prcc)),
      sprintf("%d pairs; %d not estimable", nrow(sens), sum(!is.finite(sens$prcc))))
check("the sample is large enough to separate the parameters",
      min(sens$n) - length(PARS) >= 200,
      sprintf("at least %d sets with a finite statistic and %d parameters", min(sens$n), length(PARS)))
check("at least some parameters move the calibration statistics", length(INFLUENTIAL) >= 5,
      sprintf("%d influential, %d negligible", length(INFLUENTIAL), length(NEGLIGIBLE)))

# --- identifiability: how much of each prior the region meeting R1 to R3 keeps -----------------
span_of <- function(pn) {
  v <- feas[[pn]]; v <- v[is.finite(v)]
  # the population parameters span the grid they are chosen from
  pr <- if (pn %in% names(QSP2_PRIOR)) QSP2_PRIOR[[pn]] else c(range(QSP2_TRI_GRID[[pn]]), "unif")
  lo0 <- as.numeric(pr[1]); hi0 <- as.numeric(pr[2])
  lo <- quantile(v, 0.05, names = FALSE); hi <- quantile(v, 0.95, names = FALSE)
  span <- if (pr[3] == "log") (log(hi) - log(lo)) / (log(hi0) - log(lo0)) else (hi - lo) / (hi0 - lo0)
  data.table(parameter = pn, prior_lo = lo0, prior_hi = hi0, lo = lo, hi = hi, span_frac = span,
             identified = is.finite(span) && span < 0.6)
}
prof <- rbindlist(lapply(PARS, span_of))
fwrite(prof, file.path(OUT_DIR, "B07_profiles.csv"))
cat("\nhow much of each prior survives where R1 to R3 hold (smaller means better pinned down):\n")
print(prof[order(span_frac)][, .(parameter, lo = signif(lo, 3), hi = signif(hi, 3),
                                 span = round(span_frac, 3), identified)])

check("every parameter has a reported range over the region meeting R1 to R3",
      all(is.finite(prof$span_frac)), sprintf("%d of %d", sum(is.finite(prof$span_frac)), nrow(prof)))
check("the report separates parameters the region pins down from those it does not",
      nrow(prof) > 0, sprintf("%d pinned down, %d not", sum(prof$identified), sum(!prof$identified)))

save_json(list(block = "B07",
               evaluated = nrow(all6), feasible = nrow(feas), influential = INFLUENTIAL,
               negligible = NEGLIGIBLE, identified = prof[identified == TRUE]$parameter,
               not_identified = prof[identified == FALSE]$parameter,
               fixed = character(0),
               note = "no fitted parameter has a literature value, so none is fixed; all keep their prior range in B8"),
          "B07_sensitivity.json")
