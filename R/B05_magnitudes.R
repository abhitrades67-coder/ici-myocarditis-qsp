# ---------------------------------------------------------------------------
# B05_magnitudes.R  -  block B5: the magnitude audit, before any fitting.
#
# Every observable the model produces is compared with an independent absolute
# number.  A unit error in an observable can survive fitting, because the fitted
# parameters absorb it; this block looks for one before anything is fitted.
#
# Nothing here is fitted to a calibration target.  What this block does set is the
# troponin release coefficient theta, which the specification says is chosen so that
# the model reproduces the measured peak per gram; that is a unit conversion of the
# observation model, not a fit to an endpoint.
#
# Gate: every audited magnitude within a factor of 3 of its anchor, or stopped.
#
# Outputs: output/B05_magnitudes.csv, output/B05_nominal_pars.csv, output/B05_magnitudes.json
# ---------------------------------------------------------------------------

source(file.path(QSP2_DIR, "R", "qsp2_model.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "run_tests.R"), local = TRUE)
build_qsp2_dll(); load_qsp2_dll()

cat("B5  magnitude audit\n\n")
FACTOR <- 3
within <- function(x, anchor, f = FACTOR) is.finite(x) && x > 0 && x / anchor <= f && anchor / x <= f
audit <- data.table(quantity = character(), model = numeric(), anchor = numeric(),
                    ratio = numeric(), unit = character(), source = character())
note <- function(q, m, a, u, s) audit <<- rbind(audit, data.table(
  quantity = q, model = m, anchor = a, ratio = m / a, unit = u, source = s))

tr <- qsp2_run_tests()
check("the model still passes its regression tests", isTRUE(tr$ok),
      sprintf("%d tests, %d failed, %d errors", tr$n, tr$failed, tr$errors))

pars <- qsp2_default_pars(assay = "hs_cTnI")   # Mewton measured troponin I

# =====================================================================================
cat("1. Troponin: the release coefficient and what it implies\n")
# The troponin submodule is integrated on its own against a prescribed injury, so that
# the audit is of the observation model and not of the immune model.
#
#   dQs/dt = (1 - fc) theta k(t) Mlv - k_rel Qs
#   dTn/dt = [fc theta k(t) Mlv + k_rel Qs]/Vd - k_el (Tn - Tn0)
#
# The prescribed injury destroys a given fraction of the ventricle over a given number
# of days, which is how acute myocarditis presents in the source cohort.
tn_peak <- function(theta, frac = 0.13, days = 4, p = pars) {
  k <- frac / days
  f <- function(t, y, .) {
    kap <- if (t <= days) k else 0
    rel <- theta * kap * p[["Mlv"]]
    list(c((1 - p[["fc"]]) * rel - p[["k_rel"]] * y[1],
           (p[["fc"]] * rel + p[["k_rel"]] * y[1]) / p[["Vd"]] - p[["k_el"]] * (y[2] - p[["Tn0"]])))
  }
  o <- deSolve::lsoda(c(Qs = 0, Tn = p[["Tn0"]]), seq(0, 60, by = 0.02), f, NULL,
                      rtol = 1e-10, atol = 1e-12)
  max(o[, "Tn"]) - p[["Tn0"]]
}

# solve for theta so that a Mewton-like injury peaks at 8.5 ug/L, which is 8500 ng/L
target_peak <- 8500
theta <- uniroot(function(th) tn_peak(th) - target_peak,
                 interval = c(1e3, 1e9), tol = 1e-6)$root
pars[["theta"]] <- theta

peak_per_gram <- tn_peak(theta) / (0.13 * pars[["Mlv"]])
note("peak troponin per gram of injury", peak_per_gram, 527, "ng/L per gram",
     "Mewton 2015: 8.5 ug/L peak with 13% of a 124 g ventricle")
check("the troponin release coefficient reproduces the measured peak per gram",
      within(peak_per_gram, 527),
      sprintf("theta = %.3g ng per gram; model gives %.0f ng/L per gram against 527",
              theta, peak_per_gram))

# how the coefficient depends on how fast the injury happens: reported, not hidden
th_fast <- uniroot(function(th) tn_peak(th, days = 1) - target_peak, c(1e3, 1e9))$root
th_slow <- uniroot(function(th) tn_peak(th, days = 14) - target_peak, c(1e3, 1e9))$root
check("the release coefficient is reported with its dependence on injury duration",
      is.finite(th_fast) && is.finite(th_slow),
      sprintf("injury over 1 day: %.3g; over 4 days: %.3g; over 14 days: %.3g ng per gram; a %.1f-fold span",
              th_fast, theta, th_slow, th_slow / th_fast))

# the alternative coefficient from a preprint, recorded for comparison
th_pre <- uniroot(function(th) tn_peak(th) / (0.13 * pars[["Mlv"]]) - 909, c(1e3, 1e10))$root
check("the alternative coefficient from the preprint is recorded",
      is.finite(th_pre), sprintf("909 ng/L per gram corresponds to theta = %.3g ng per gram", th_pre))

# smallest injury that a measurement can see
lim_I <- QSP2_TROPONIN$hs_cTnI$limit
g_detect <- uniroot(function(g) tn_peak(theta, frac = g / pars[["Mlv"]]) -
                      (lim_I - pars[["Tn0"]]), c(1e-6, 5))$root
note("myocardium destroyed to reach the assay limit", g_detect, 0.0497, "g",
     "arithmetic from the specification, Section 9")
check("the smallest visible injury is a small fraction of the ventricle",
      within(g_detect, 0.0497), sprintf("%.4f g, which is %.4f%% of the ventricle",
                                        g_detect, 100 * g_detect / pars[["Mlv"]]))

# =====================================================================================
cat("\n2. Pharmacokinetic exposure against published anchors\n")
BW <- 80
pp <- qsp2_default_pars(BW = BW, assay = "hs_cTnI"); pp[["theta"]] <- theta
ipi <- qsp2_simulate(pp, qsp2_regimen("ipi_alone", BW = BW), seq(0, 90, by = 0.25))
cmin1 <- ipi$Ci[which.min(abs(ipi$t - 21))]
note("ipilimumab trough after the first 3 mg/kg dose", cmin1, 13, "ug/mL",
     "Feng 2014: about 90% of patients above 6 ug/mL, few above 20")
check("the ipilimumab trough after one dose sits in the published range",
      cmin1 > 6 && cmin1 < 20,
      sprintf("%.1f ug/mL; the source says about 90%% exceed 6 and about 2%% exceed 20", cmin1))

nv <- qsp2_simulate(pp, qsp2_regimen("nivo_q2w"), seq(0, 200, by = 0.25))
cmax_n <- max(nv$Cn); cmin_n <- nv$Cn[which.min(abs(nv$t - 196))]
note("nivolumab peak concentration at steady state", cmax_n, 100, "ug/mL",
     "240 mg into a central volume of 3.63 L, with accumulation")
check("nivolumab concentrations are of the expected order",
      cmax_n > 20 && cmax_n < 300 && cmin_n > 5,
      sprintf("peak %.1f ug/mL, trough %.1f ug/mL", cmax_n, cmin_n))
check("PD-1 is close to saturated at label doses",
      max(nv$Opd1) > 0.95 && nv$Opd1[which.min(abs(nv$t - 196))] > 0.9,
      sprintf("occupancy peaks at %.3f and troughs at %.3f",
              max(nv$Opd1), nv$Opd1[which.min(abs(nv$t - 196))]))

# =====================================================================================
cat("\n3. Tissue loss and ejection fraction\n")
lvef <- function(v, d, p) qsp2_lvef(v, d, p)
ef_10 <- lvef(0.90, 0.10, pars); ef_25 <- lvef(0.75, 0.25, pars); ef_0 <- lvef(1, 0, pars)
note("ejection fraction with a tenth of the ventricle damaged", ef_10, 50, "%",
     "the registry's cut-off for reduced function")
check("the function map gives a sensible ejection fraction at baseline and with damage",
      abs(ef_0 - 60) < 1e-9 && ef_10 > 45 && ef_10 < 60 && ef_25 < ef_10,
      sprintf("baseline %.1f%%, a tenth damaged %.1f%%, a quarter damaged %.1f%%",
              ef_0, ef_10, ef_25))

# what loss is needed to fall below the registry cut-off
f_below <- uniroot(function(x) lvef(1 - x, x * 0.5, pars) - 50, c(1e-4, 0.9))$root
note("loss needed to bring ejection fraction below 50 percent", f_below, 0.20, "fraction",
     "judgement: a fifth of the ventricle is the expected order")
check("a plausible amount of damage is needed to reach the severe tier",
      within(f_below, 0.20), sprintf("%.1f%% of the ventricle", 100 * f_below))

# =====================================================================================
cat("\n4. Immune magnitudes\n")
# The tier T2 threshold in cells, from the specification's conversion
T2_cells <- QSP2_T2$cells
note("T cells in the ventricle at the biopsy threshold", T2_cells, 7.5e7, "cells",
     "Caforio 2013 threshold with the Abercrombie correction")
check("the biopsy threshold in cells matches the specification", within(T2_cells, 7.5e7, 1.05),
      sprintf("%.2e cells", T2_cells))

# can the model reach it, and is the resting heart far below it?
pk <- pars; pk[["Sus"]] <- 20        # a susceptible patient, to test what is reachable
sim <- qsp2_simulate(pk, qsp2_regimen("combo_nivo1_ipi3"), seq(0, 365, by = 1))
note("resting cardiac T cells as a share of the threshold",
     pars[["Tn0"]] * 0 + qsp2_init_state(pars)[["Rm"]] / T2_cells, 0.002, "fraction",
     "a healthy heart must sit well below the threshold")
check("a resting heart is far below the biopsy threshold",
      qsp2_init_state(pars)[["Rm"]] < 0.05 * T2_cells,
      sprintf("%.2e resting cells against a threshold of %.2e", qsp2_init_state(pars)[["Rm"]], T2_cells))
check("the model can reach the biopsy threshold during a response",
      max(sim$Tcells_heart) > T2_cells,
      sprintf("peak %.2e cells against a threshold of %.2e", max(sim$Tcells_heart), T2_cells))

# The tier is f_dx times this, with f_dx fitted from 0.03 upwards.  The
# bottom of that prior has to leave a heart with no disease clearly below the tier,
# otherwise the model would call a healthy resident-memory pool a case of myocarditis.
F_DX_MIN <- 0.03
check("even at the lowest myocarditis threshold the prior allows, a resting heart is not a case",
      qsp2_init_state(pars)[["Rm"]] < 0.05 * F_DX_MIN * T2_cells,
      sprintf("%.2e resting cells against %.2e, the tier at f_dx = %.2f",
              qsp2_init_state(pars)[["Rm"]], F_DX_MIN * T2_cells, F_DX_MIN))

# =====================================================================================
cat("\n5. Are the nominal outcomes physiological?\n")
# The default parameters are mid-prior guesses, not a calibration.  What this checks is
# that the nominal setting does not destroy the heart, which would mean the killing
# term is out by orders of magnitude before fitting ever starts.
worst_loss <- 1 - min(sim$v)
check("the nominal parameters do not destroy the ventricle",
      worst_loss < 0.5,
      sprintf("%.1f%% of the ventricle lost at nominal parameters with susceptibility 20",
              100 * worst_loss))

fwrite(audit, file.path(OUT_DIR, "B05_magnitudes.csv"))
np <- data.table(parameter = names(pars), value = as.numeric(pars))
fwrite(np, file.path(OUT_DIR, "B05_nominal_pars.csv"))
cat(sprintf("\n%d magnitudes audited, recorded in output/B05_magnitudes.csv\n", nrow(audit)))
print(audit[, .(quantity, model = signif(model, 4), anchor = signif(anchor, 4),
                ratio = round(ratio, 3))])

save_json(list(block = "B05",
               factor_allowed = FACTOR, theta = theta,
               theta_injury_1day = th_fast, theta_injury_14day = th_slow,
               theta_preprint_alternative = th_pre,
               peak_ng_per_L_per_gram = peak_per_gram,
               grams_to_reach_assay_limit = g_detect,
               T2_threshold_cells = T2_cells,
               worst_nominal_tissue_loss = worst_loss,
               audited = nrow(audit)), "B05_magnitudes.json")
