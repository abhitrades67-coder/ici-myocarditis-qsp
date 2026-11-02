# ---------------------------------------------------------------------------
# B04_implementation.R  -  block B4: the model is implemented and verified.
#
# The model is src/qsp2_rhs.c with the R layer in R/qsp2_model.R.  This block
# does not fit anything and does not look at any target.  It asks one question:
# does the code do what docs/model_specification.md says it does?
#
# The checks fall into nine groups:
#   1. structure           states, parameters, and the two implementations agree
#   2. conservation        tissue adds up, scar never falls, nothing goes negative
#   3. analytic limits     no drug, no killing, no troponin
#   4. pharmacokinetics    half-lives and accumulation against the published values
#   5. numerics            solver converged, dose delivered in full, units consistent
#   6. structure of the answer   the latent phase is real, not imposed
#   7. observation model   assay, tier and threshold arguments do what they say
#   8. population engine   tempo, treatment, noise, baseline, background, severity,
#                          and the engine against a Monte Carlo of single patients
#   9. regression tests    the testthat suite in tests/
#
# Outputs: output/B04_verification.csv, output/B04_engine_vs_mc.csv,
#          output/B04_implementation.json
# ---------------------------------------------------------------------------

source(file.path(QSP2_DIR, "R", "qsp2_model.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "qsp2_rhs_reference.R"), local = TRUE)
build_qsp2_dll(); load_qsp2_dll()

cat("B4  implementation and verification\n\n")
approx <- function(a, b, tol) all(abs(a - b) <= tol * abs(b))

# =====================================================================================
cat("1. Structure\n")
pars <- qsp2_default_pars()
y0   <- qsp2_init_state(pars)
NT   <- as.integer(round(pars[["NT"]]))

check("the parameter vector has the length the C code expects",
      all(QSP2_PAR_NAMES %in% names(pars)) &&
        identical(names(pars[QSP2_PAR_NAMES]), QSP2_PAR_NAMES),
      sprintf("%d parameters, of which %d go to the solver",
              length(pars), length(QSP2_PAR_NAMES)))
check("the state vector has 19 states plus the transit chain",
      length(y0) == 19L + NT, sprintf("%d states with NT = %d", length(y0), NT))
check("no parameter is missing or not finite", all(is.finite(pars)),
      paste(names(pars)[!is.finite(pars)], collapse = ", "))

# Every term the specification writes into an equation must exist somewhere in the code.
# Comparing the two right-hand sides with each other does not show that either matches the
# document, so the document is read as well.
spec_txt <- paste(readLines(file.path(DOC_DIR, "model_specification.md"), warn = FALSE,
                            encoding = "UTF-8"), collapse = "
")
code_txt <- paste(c(readLines(file.path(SRC_DIR, "qsp2_rhs.c"), warn = FALSE),
                    readLines(file.path(QSP2_DIR, "R", "qsp2_model.R"), warn = FALSE),
                    readLines(file.path(QSP2_DIR, "R", "lib_vpop.R"), warn = FALSE),
                    readLines(file.path(QSP2_DIR, "R", "lib_pop.R"), warn = FALSE)),
                  collapse = "
")
TERMS <- c(background = "background", theta = "theta", k_rel = "k_rel", fc = "fc",
           k_el = "k_el", kappa = "kappa", licensing = "k_lic", exhaustion = "eps_x",
           transit = "k_tr", occupancy = "EC50pd1", ceiling = "Emax", priming = "lam_p",
           vulnerability = "Svul", symptom = "k_sym")
present <- vapply(TERMS, function(x) grepl(x, code_txt, fixed = TRUE), logical(1))
in_spec <- vapply(TERMS, function(x) grepl(x, spec_txt, fixed = TRUE), logical(1))
missing_terms <- names(TERMS)[in_spec & !present]
check("every term the specification puts in an equation exists in the code",
      !length(missing_terms),
      if (length(missing_terms)) paste("in the specification but not in the code:",
                                       paste(missing_terms, collapse = ", "))
      else sprintf("%d terms checked", sum(in_spec)))
# The effector pool must contract whenever no receptor is blocked.  Without this the
# untreated model generates autoimmunity.
check("the effector pool contracts when no checkpoint is blocked",
      pars[["rho"]] < pars[["del_E"]] + pars[["k_tr"]],
      sprintf("rho = %.3f against del_E + k_tr = %.3f", pars[["rho"]],
              pars[["del_E"]] + pars[["k_tr"]]))
check("with no drug the priming term is exactly zero",
      qsp2_rhs_R(0, y0, pars)[[1]][7] == 0,
      "dP/dt at baseline with no drug")

# The central check of this block: two implementations, written separately from the
# same document, must agree.  A term typed wrongly in either one shows up here.
set.seed(20260927)
nrand <- 1000L
worst <- 0; worst_state <- NA_integer_
for (i in seq_len(nrand)) {
  yr <- y0
  yr[1:6]  <- runif(6, 0, 500)                       # drug amounts, mg
  yr[7:8]  <- 10^runif(2, 0, 9)                      # P and E
  yr[9:(8 + NT)] <- 10^runif(NT, 0, 8)               # transit chain
  j <- 8L + NT
  yr[j + 1] <- 10^runif(1, 0, 8)                     # Eh
  yr[j + 2] <- 10^runif(1, 0, 6)                     # Rm
  yr[j + 3] <- 10^runif(1, 0, 6)                     # L
  yr[j + 4] <- 10^runif(1, 0, 7)                     # Mq
  yr[j + 5] <- 10^runif(1, 0, 6)                     # Gr
  yr[j + 6] <- runif(1, 0, 50)                       # X
  vv <- runif(1, 0, 1); dd <- runif(1, 0, 1 - vv)
  yr[j + 7] <- vv; yr[j + 8] <- dd; yr[j + 9] <- 1 - vv - dd
  yr[j + 10] <- 10^runif(1, 0, 8)                    # Qs
  yr[j + 11] <- 10^runif(1, 0, 4)                    # Tn
  pr <- pars
  pr[["infn"]] <- runif(1, 0, 2000); pr[["infi"]] <- runif(1, 0, 2000)
  pr[["infp"]] <- runif(1, 0, 2000)
  a <- qsp2_rhs_C(runif(1, 0, 365), yr, pr)$dy
  b <- qsp2_rhs_R(0, yr, pr)[[1]]
  sc <- pmax(abs(a), abs(b), 1e-30)
  rel <- max(abs(a - b) / sc)
  if (rel > worst) { worst <- rel; worst_state <- i }
}
check("the compiled and the independent R right-hand sides agree to 1e-10",
      worst < 1e-10,
      sprintf("worst relative difference %.3e over %d random states (state %d)",
              worst, nrand, worst_state))

# =====================================================================================
cat("\n2. Conservation\n")
reg  <- qsp2_regimen("combo_nivo1_ipi3")
tt   <- seq(0, 365, by = 0.5)
pk   <- pars; pk[["k_kill"]] <- 0.05; pk[["Sus"]] <- 3
run  <- qsp2_simulate(pk, reg, tt)

tot <- run$v + run$d + run$s
check("viable, damaged and scar always add to one", max(abs(tot - 1)) < 1e-9,
      sprintf("worst departure %.2e", max(abs(tot - 1))))
check("scar never decreases", all(diff(run$s) >= -1e-12),
      sprintf("smallest change in scar %.2e", min(diff(run$s))))
statecols <- intersect(names(run), names(y0))
mins <- vapply(statecols, function(cn) min(run[[cn]]), numeric(1))
check("no state goes negative", all(mins >= -1e-9),
      sprintf("smallest value %.2e in %s", min(mins), statecols[which.min(mins)]))
check("receptor occupancy stays between zero and one",
      min(run$Opd1, run$Octla4) >= 0 && max(run$Opd1, run$Octla4) <= 1,
      sprintf("PD-1 in [%.3f, %.3f], CTLA-4 in [%.3f, %.3f]",
              min(run$Opd1), max(run$Opd1), min(run$Octla4), max(run$Octla4)))

# =====================================================================================
cat("\n3. Analytic limits\n")
none <- qsp2_simulate(pars, qsp2_regimen("none"), seq(0, 730, by = 1))
moved <- vapply(statecols, function(cn) max(abs(none[[cn]] - none[[cn]][1])), numeric(1))
check("with no drug nothing moves over two years", max(moved) < 1e-6,
      sprintf("largest movement %.2e in %s", max(moved), statecols[which.max(moved)]))

p0 <- pars; p0[["k_kill"]] <- 0
z  <- qsp2_simulate(p0, reg, tt)
check("with no killing there is no injury", max(abs(z$v - 1)) < 1e-12,
      sprintf("viable tissue departs from one by at most %.2e", max(abs(z$v - 1))))
check("with no killing troponin stays at baseline",
      max(abs(z$Tn - pars[["Tn0"]])) < 1e-9,
      sprintf("troponin departs from baseline by at most %.2e", max(abs(z$Tn - pars[["Tn0"]]))))
check("with no drug there is no receptor occupancy",
      max(none$Opd1, none$Octla4) < 1e-12, "occupancy is zero throughout")

# troponin alone must decay at the rate its half-life implies
yt <- y0; yt[["Tn"]] <- pars[["Tn0"]] + 1000
dec <- qsp2_simulate(p0, qsp2_regimen("none"), c(0, 0.25), y0 = yt)
pred <- pars[["Tn0"]] + 1000 * exp(-pars[["k_el"]] * 0.25)
check("troponin decays at the rate its measured half-life implies",
      approx(dec$Tn[2], pred, 1e-6),
      sprintf("after 6 hours: model %.4f, expected %.4f", dec$Tn[2], pred))

# =====================================================================================
cat("\n4. Pharmacokinetics against the published values\n")
# For a two-compartment model the two half-lives follow from the parameters, so they
# are computed here in closed form and compared with the papers.
hl <- function(CL, Vc, Q, Vp) {
  k10 <- CL / Vc; k12 <- Q / Vc; k21 <- Q / Vp
  s <- k10 + k12 + k21; r <- sqrt(s^2 - 4 * k10 * k21)
  c(alpha = log(2) / (0.5 * (s + r)), beta = log(2) / (0.5 * (s - r)))
}
h_n <- hl(pars[["CLn"]], pars[["Vcn"]], pars[["Qn"]], pars[["Vpn"]])
h_i <- hl(pars[["CLi"]], pars[["Vci"]], pars[["Qi"]], pars[["Vpi"]])
h_p <- hl(pars[["CLp"]], pars[["Vcp"]], pars[["Qp"]], pars[["Vpp"]])

check("ipilimumab terminal half-life matches the published 14.7 days",
      approx(h_i[["beta"]], 14.7, 0.05), sprintf("model %.2f days", h_i[["beta"]]))
check("pembrolizumab terminal half-life matches the published 27.3 days",
      approx(h_p[["beta"]], 27.3, 0.08), sprintf("model %.2f days", h_p[["beta"]]))
check("nivolumab distribution half-life matches the published 32.5 hours",
      approx(h_n[["alpha"]] * 24, 32.5, 0.05), sprintf("model %.1f hours", h_n[["alpha"]] * 24))
check("ipilimumab distribution half-life matches the published 27.4 hours",
      approx(h_i[["alpha"]] * 24, 27.4, 0.05), sprintf("model %.1f hours", h_i[["alpha"]] * 24))

# Nivolumab is the one case where the specification's simplification bites: the
# published 25 days is measured at steady state, where clearance has fallen by about
# 24.5%.  With constant clearance the model must sit below it, and with clearance
# reduced by 24.5% it must sit above it.  The published value has to lie between.
h_n_ss <- hl(pars[["CLn"]] * (1 - 0.245), pars[["Vcn"]], pars[["Qn"]], pars[["Vpn"]])
check("the published nivolumab half-life lies between the model's two clearance cases",
      h_n[["beta"]] < 25 && h_n_ss[["beta"]] > 25,
      sprintf("constant clearance %.1f days, clearance reduced by 24.5%% %.1f days, published 25 days",
              h_n[["beta"]], h_n_ss[["beta"]]))

# accumulation of pembrolizumab with three-weekly dosing
acc <- qsp2_simulate(pars, qsp2_regimen("pembro_q3w"), seq(0, 400, by = 0.5))
tr1 <- acc$Cp[which.min(abs(acc$t - 21))]
trss <- acc$Cp[which.min(abs(acc$t - 378))]
check("pembrolizumab accumulates about 2.2-fold with three-weekly dosing",
      approx(trss / tr1, 2.2, 0.20),
      sprintf("trough ratio %.2f (published about 2.2-fold)", trss / tr1))

# the unit conversions written into the model
check("the pharmacokinetic unit conversions are exact",
      approx(QSP2_PK$nivolumab$CL, 9.4 * 24 / 1000, 1e-12) &&
        approx(QSP2_PK$nivolumab$Q, 32.1 * 24 / 1000, 1e-12) &&
        approx(QSP2_PK$ipilimumab$CL, 0.0150 * 24, 1e-12) &&
        approx(QSP2_PK$ipilimumab$Q, 0.0411 * 24, 1e-12),
      sprintf("nivolumab CL %.4f L/day, ipilimumab CL %.4f L/day",
              QSP2_PK$nivolumab$CL, QSP2_PK$ipilimumab$CL))

# =====================================================================================
cat("\n5. Numerics\n")
# Every milligram must arrive.  With no elimination, the amount in the body at the end
# has to equal the amount given.
pz <- pars
for (nm in c("CLn", "CLi", "CLp")) pz[[nm]] <- 0
rz <- qsp2_regimen("combo_nivo1_ipi3")
sz <- qsp2_simulate(pz, rz, c(0, 365))
given_n <- sum(rz[drug == "nivolumab"]$dose)
given_i <- sum(rz[drug == "ipilimumab"]$dose)
in_body_n <- sz$A1n[2] + sz$A2n[2]; in_body_i <- sz$A1i[2] + sz$A2i[2]
check("every dose is delivered in full across the dose boundaries",
      approx(in_body_n, given_n, 1e-6) && approx(in_body_i, given_i, 1e-6),
      sprintf("nivolumab given %.1f mg, in body %.4f mg; ipilimumab given %.1f mg, in body %.4f mg",
              given_n, in_body_n, given_i, in_body_i))

# solver tolerance
loose <- qsp2_simulate(pk, reg, tt, rtol = 1e-8,  atol = 1e-10)
tight <- qsp2_simulate(pk, reg, tt, rtol = 1e-12, atol = 1e-14)
d_tol <- max(abs(loose$Tn - tight$Tn) / pmax(abs(tight$Tn), 1e-10))
check("the answer does not depend on the solver tolerance", d_tol < 1e-5,
      sprintf("largest relative difference in troponin %.2e between 1e-8 and 1e-12", d_tol))

# output grid
coarse <- qsp2_simulate(pk, reg, seq(0, 365, by = 5))
fine   <- run[t %in% seq(0, 365, by = 5)]
d_grid <- max(abs(coarse$Tn - fine$Tn) / pmax(abs(fine$Tn), 1e-10))
check("the answer does not depend on the output grid", d_grid < 1e-6,
      sprintf("largest relative difference %.2e between a 5 day and a 0.5 day grid", d_grid))

# time in hours instead of days: a dimensional check
rate_pars <- c("CLn", "Qn", "CLi", "Qi", "CLp", "Qp", "lam_p", "k_act", "del_P", "rho",
               "del_E", "k_tr", "pi_h", "del_h", "eps_x", "k_lic", "del_L", "k_m", "del_m",
               "lam_g", "k_g", "del_g", "del_x", "k_kill", "r_rep", "r_scar", "k_rel", "k_el")
ph <- pk; ph[rate_pars] <- ph[rate_pars] / 24
rh <- copy(reg); rh[, `:=`(day = day * 24, dur = dur * 24)]
hrs <- qsp2_simulate(ph, rh, tt * 24)
d_units <- max(abs(hrs$Tn - run$Tn) / pmax(abs(run$Tn), 1e-10))
check("running the model in hours instead of days gives the same answer", d_units < 1e-6,
      sprintf("largest relative difference in troponin %.2e", d_units))

# =====================================================================================
cat("\n6. The latent phase is a property of the model, not an imposed rule\n")
# Nothing in the code forbids early killing.  It cannot happen because resident memory
# cells are unlicensed until effectors arrive, and effectors have to be made and travel.
first_pos <- function(x, tv, eps = 1e-12) if (any(x > eps)) tv[which(x > eps)[1]] else Inf
t_L  <- first_pos(run$L, run$t)
t_Eh <- first_pos(run$Eh, run$t)
t_inj <- first_pos(1 - run$v, run$t, 1e-9)
check("licensed resident memory cells appear only after effectors reach the heart",
      t_L >= t_Eh, sprintf("effectors at day %.1f, licensed cells at day %.1f", t_Eh, t_L))
# How long the latent phase is depends on the parameters and is requirement R1, which
# block B6 tests across the whole prior.  What B4 checks is that the structure forces a
# latent phase at all: injury cannot start before effectors have arrived in the heart.
check("injury cannot begin before effectors reach the heart", t_inj > t_Eh,
      sprintf("effectors at day %.1f, first measurable injury at day %.1f (the length of the latent phase is tested in B6)",
              t_Eh, t_inj))
check("at the first dose nothing in the heart is able to kill",
      run$kappa[run$t == 0] == 0, "the killing rate at time zero is exactly zero")

# Requirement R5's precondition: with ipilimumab alone the effector pool still grows.
ipi <- qsp2_simulate(pk, qsp2_regimen("ipi_alone"), seq(0, 365, by = 1))
check("with ipilimumab alone the effector pool still expands",
      max(ipi$E) > 10 * ipi$E[1] + 1,
      sprintf("effector pool reaches %.3e from %.3e", max(ipi$E), ipi$E[1]))
check("with ipilimumab alone there is CTLA-4 blockade and no PD-1 blockade",
      max(ipi$Octla4) > 0.5 && max(ipi$Opd1) < 1e-12,
      sprintf("CTLA-4 occupancy peaks at %.3f, PD-1 at %.2e", max(ipi$Octla4), max(ipi$Opd1)))

# dosing calendars
r1 <- qsp2_regimen("nivo_q2w"); r2 <- qsp2_regimen("combo_nivo1_ipi3")
check("the dosing calendars match the prescribing information",
      all(diff(r1$day) == 14) && nrow(r2[drug == "ipilimumab"]) == 4 &&
        all(r2[drug == "ipilimumab"]$day == c(0, 21, 42, 63)),
      sprintf("nivolumab every 14 days; ipilimumab on days %s",
              paste(r2[drug == "ipilimumab"]$day, collapse = ", ")))

# =====================================================================================
cat("
7. The observation model says what it is asked to say
")
# Each check tests by effect that an argument or a constant does what the caller believes
# it does.
source(file.path(QSP2_DIR, "R", "lib_vpop.R"), local = TRUE)
pv <- qsp2_default_pars(assay = "hs_cTnI")
pv[["Sus"]] <- 20
regv <- qsp2_regimen("combo_nivo1_ipi3"); ttv <- seq(0, 200, by = 1)

# (a) the assay argument must change the troponin kinetics
si <- qsp2_simulate(pv, regv, ttv, assay = "hs_cTnI")
st <- qsp2_simulate(pv, regv, ttv, assay = "hs_cTnT")
check("asking for a different assay changes the troponin trajectory",
      max(abs(si$Tn - st$Tn)) > 1e-6 * max(si$Tn),
      sprintf("peak %.4g ng/L with troponin I kinetics against %.4g with troponin T",
              max(si$Tn), max(st$Tn)))
check("the faster-clearing assay gives the lower peak",
      max(st$Tn) < max(si$Tn),
      sprintf("half-life %.0f min for T against %.0f min for I",
              log(2) / QSP2_TROPONIN$hs_cTnT$k_el * 1440,
              log(2) / QSP2_TROPONIN$hs_cTnI$k_el * 1440))

# (b) the myocarditis tier has to move when f_dx moves
pl <- pv; pl[["f_dx"]] <- 0.1
check("the myocarditis threshold scales with the fitted factor",
      abs(qsp2_t2_cells(pl) - 0.1 * QSP2_T2$cells) < 1e-6 * QSP2_T2$cells &&
        abs(qsp2_t2_cells(pv) - QSP2_T2$cells) < 1e-6 * QSP2_T2$cells,
      sprintf("%.3e cells at f_dx = 1, %.3e at f_dx = 0.1",
              qsp2_t2_cells(pv), qsp2_t2_cells(pl)))
cv1 <- qsp2_response_curve(pv, regv, times = ttv, schedules = "S2",
                           sgrid = 10^seq(-2, 3, by = 0.5))
cv2 <- qsp2_response_curve(pl, regv, times = ttv, schedules = "S2",
                           sgrid = 10^seq(-2, 3, by = 0.5))
check("a lower myocarditis threshold makes the tier easier to reach, never harder",
      sum(cv2$T2) >= sum(cv1$T2),
      sprintf("%d of %d susceptibilities reach it at f_dx = 1, %d at f_dx = 0.1",
              sum(cv1$T2), nrow(cv1), sum(cv2$T2)))

# (c) a cohort has to be scored at its own threshold
COHV <- list(low  = qsp2_cohort(seq(0, 200, by = 14), "hs_cTnI", 14),
             high = qsp2_cohort(seq(0, 200, by = 14), "hs_cTnI", 200))
cvc <- qsp2_response_curve(pv, regv, times = ttv, schedules = "S2",
                           sgrid = 10^seq(-2, 3, by = 0.5), cohorts = COHV)
check("a cohort with a higher threshold never detects more than one with a lower",
      all(cvc$det_high <= cvc$det_low),
      sprintf("%d of %d susceptibilities detected at 14 ng/L, %d at 200 ng/L",
              sum(cvc$det_low), nrow(cvc), sum(cvc$det_high)))
check("the two thresholds do not give the same answer, so the limit is being read",
      sum(cvc$det_low) > sum(cvc$det_high),
      sprintf("%d against %d", sum(cvc$det_low), sum(cvc$det_high)))

# =====================================================================================
cat("\n8. The population engine\n")
# Each check tests a mechanism by what it does, not by whether a name exists.  The engine is
# R/lib_pop.R.  R/lib_mc.R is a Monte Carlo of individual patients written separately from
# it, which follows the treated course and the background episodes explicitly.
source(file.path(QSP2_DIR, "R", "lib_pop.R"), local = TRUE)
source(file.path(QSP2_DIR, "R", "lib_mc.R"), local = TRUE)
# A hand-set parameter vector, not a calibration: chosen so that the population has
# myocarditis, severe cases and cases that screening prevents, so the comparisons can fail.
PE <- qsp2_default_pars(assay = "hs_cTnI")
for (nm in names(m8 <- list(f_dx = 0.1, k_kill = 3e-2, L_tx = 10, k_sym = 20, s_ksym = 0.8,
                            eta = 0.1, pexp = 2.5, s_b = 1, lam_bg = 0.05, m_bg = 40)))
  PE[[nm]] <- m8[[nm]]
PP8 <- qsp2_pop_pp(PE)

# (a) tempo
t1_at <- function(tau) {
  q <- qsp2_apply_tempo(PE, tau); q[["Sus"]] <- 3
  r <- qsp2_simulate(q, qsp2_regimen("nivo_q2w"), 0:365, assay = "hs_cTnI")
  w <- which(r$Tn > QSP2_REF_URL); if (length(w)) w[1] - 1 else Inf
}
tt8 <- vapply(c(0.5, 1, 2), t1_at, 0)
check("a slower tempo gives a later first troponin rise, a faster one an earlier",
      tt8[1] > tt8[2] && tt8[2] > tt8[3],
      sprintf("first day above the limit %s at tempo 0.5, 1 and 2", paste(tt8, collapse = ", ")))
pq <- qsp2_apply_tempo(PE, 2)
chg <- names(PE)[abs(pq - PE) > 1e-12 * abs(PE)]
check("tempo multiplies exactly the rates of the clone's response, and nothing else",
      setequal(chg, QSP2_TAU_RATES) &&
        all(abs(pq[QSP2_TAU_RATES] / PE[QSP2_TAU_RATES] - 2) < 1e-12) &&
        # the myocardium's own rates, the killing constant and every capacity
        # keep their values, because they describe tissue rather than the clone
        !any(c("r_rep", "r_scar", "k_kill", "K_lic", "Kg", "Mmax", "Kh") %in% QSP2_TAU_RATES),
      sprintf("%d rates scale: %s", length(QSP2_TAU_RATES), paste(chg, collapse = ", ")))
# tempo leaves the resting state untouched because the quantities that set it are ratios of
# two rates that carry the same factor
y_ref <- qsp2_init_state(PE); y_tau <- qsp2_init_state(qsp2_apply_tempo(PE, 3))
check("a patient's tempo does not move the resting state",
      max(abs(y_ref - y_tau)) < 1e-9,
      sprintf("largest difference %.1e across %d states at a tempo of 3",
              max(abs(y_ref - y_tau)), length(y_ref)))
# The end cells of qsp2_pop_cellmass run to minus and plus infinity, so a grid narrower than
# the distribution it carries piles both tails onto its two edge nodes instead of representing
# them.  The grid must therefore be wide enough for the widest sd_tau the prior allows, and the
# check exists so that widening the prior without widening the grid cannot pass unnoticed.
# the widest spread the search can choose, read from the grid the calibration actually uses
# rather than from a figure repeated here
source(file.path(QSP2_DIR, "R", "lib_cal.R"), local = TRUE)
sd_hi <- max(QSP2_TRI_GRID$sd_tau)
w_hi  <- qsp2_pop_cellmass(QSP2_TAU_GRID, 0, sd_hi)
edge_hi <- w_hi[1] + w_hi[length(w_hi)]
check("the tempo grid can hold the widest tempo spread its prior allows",
      edge_hi < 0.01,
      sprintf("at sd_tau %.2f the two edge nodes hold %.2f%% of the population; the grid spans %.0f-fold over %d nodes",
              sd_hi, 100 * edge_hi, max(QSP2_TAU_GRID) / min(QSP2_TAU_GRID), length(QSP2_TAU_GRID)))

# (b) troponin T from troponin I, against a full troponin T solve
qss <- NULL
for (S in c(0.05, 3, 30)) for (tau in c(0.5, 2)) {
  q <- qsp2_apply_tempo(PE, tau); q[["Sus"]] <- S
  ri <- qsp2_simulate(q, qsp2_regimen("combo_nivo1_ipi3"), 0:365, assay = "hs_cTnI")
  rt <- qsp2_simulate(q, qsp2_regimen("combo_nivo1_ipi3"), 0:365, assay = "hs_cTnT")
  xi <- ri$Tn - q[["Tn0"]]; xt <- rt$Tn - q[["Tn0"]]
  bT <- q[["Tn0"]] * QSP2_TNT_URL / QSP2_REF_URL
  qss <- rbind(qss, c(err = max(abs(xt - QSP2_TNT_RATIO * xi)) / max(xt),
                      lag = abs(which(xt + bT > QSP2_TNT_URL)[1] - which(QSP2_TNT_RATIO * xi + bT > QSP2_TNT_URL)[1])))
}
check("troponin T from the clearance ratio agrees with a full troponin T solve",
      max(qss[, "err"]) < 0.03 && max(qss[, "lag"]) <= 1,
      sprintf("largest error %.2f%% of the peak; first day above 14 ng/L differs by at most %d day",
              100 * max(qss[, "err"]), max(qss[, "lag"])))

# (c) treatment, measured within the follow-up: the engine's nadir after a halt, against the
# same patient's course simulated explicitly, and an earlier halt, or any treatment, never worse.
# Before the halt the course carries the inflammatory depression; from the halt the
# immunosuppression that stops killing also resolves the infiltrate, so the rest is tissue alone.
q <- PE; q[["Sus"]] <- 3
reg_h <- qsp2_regimen("combo_nivo1_ipi3")
tr1 <- qsp2_pop_traj(PE, reg_h, sgrid = 3, taugrid = 1)          # the engine, one patient
gp1 <- qsp2_pop_prep(tr1, PP8, H = 365)$gpost[1, ]
G_dep <- tr1$G[1, ]                                               # with the depression
halt <- t(vapply(c(20, 40, 60, 90, 200, 350), function(th) {
  rh <- qsp2_simulate(q, reg_h, 0:365, assay = "hs_cTnI")
  y <- unlist(rh[t == th, names(qsp2_init_state(q)), with = FALSE])
  q0 <- q; q0[["k_kill"]] <- 0
  rc <- qsp2_simulate(q0, qsp2_regimen("none"), seq(th, 365, by = 1), y0 = y, assay = "hs_cTnI")
  c(explicit = min(min(G_dep[1:(th + 1)]), min((rc$v + q[["eta"]] * rc$d)^q[["pexp"]])),
    engine = gp1[th + 1])
}, numeric(2)))
check("the engine's nadir after a halt matches the treated course simulated in full",
      max(abs(halt[, "explicit"] - halt[, "engine"])) < 1e-6,
      sprintf("largest difference %.1e over halts on days 20, 40, 60, 90, 200 and 350, up to day 365",
              max(abs(halt[, "explicit"] - halt[, "engine"]))))
check("an earlier halt never gives a lower nadir, and treatment never a lower one than none",
      all(diff(halt[, "engine"]) <= 1e-12) && all(halt[, "engine"] >= min(G_dep) - 1e-12),
      sprintf("nadir of G %s; untreated %.3f", paste(sprintf("%.3f", halt[, "engine"]), collapse = ", "),
              min(G_dep)))

# (c2) inflammation depresses function reversibly.  A gentler patient than the
# one above, whose ventricle survives, so the depression is visible rather than hidden by the floor.
PI <- qsp2_default_pars(assay = "hs_cTnI")
for (nm in names(mI <- list(f_dx = 0.1, k_kill = 3e-4, eta = 0.5, pexp = 2.5, phi_inf = 0.4)))
  PI[[nm]] <- mI[[nm]]
PI[["Sus"]] <- 3
ri <- qsp2_simulate(PI, reg_h, 0:365, assay = "hs_cTnI")
sat_i <- qsp2_inflam_sat(ri, PI)
ef_tis <- qsp2_lvef(ri$v, ri$d, PI)
ef_inf <- qsp2_lvef(ri$v, ri$d, PI, sat = sat_i)
check("with no inflammatory burden the function map is the tissue-only one",
      max(abs(qsp2_lvef(ri$v, ri$d, PI, sat = 0) - ef_tis)) < 1e-12 &&
        abs(qsp2_lvef(1, 0, PI, sat = 0) - QSP2_MYO$EF0) < 1e-12,
      "sat = 0 reproduces the tissue-only map, and an intact ventricle sits at the baseline")
check("inflammation lowers the ejection fraction without any more tissue being lost",
      min(ef_inf) < min(ef_tis) - 1 && min(ri$v) > 0.5,
      sprintf("nadir %.1f%% from tissue alone, %.1f%% with inflammation, at a loss of %.1f%% of the ventricle",
              min(ef_tis), min(ef_inf), 100 * (1 - min(ri$v))))
# Reversibility is tested on a finite course (the four induction doses) followed for two years,
# because a patient kept on the drug all year keeps the infiltrate that depresses function.
rrev <- qsp2_simulate(PI, qsp2_regimen("combo_nivo1_ipi3", horizon = 63), 0:730, assay = "hs_cTnI")
srev <- qsp2_inflam_sat(rrev, PI)
efrev <- qsp2_lvef(rrev$v, rrev$d, PI, sat = srev)
eftrev <- qsp2_lvef(rrev$v, rrev$d, PI)
jpk <- which.max(srev)
jrec <- which(srev < 0.5 * srev[jpk] & seq_along(srev) > jpk)
check("the depression is reversible: function returns as the infiltrate contracts",
      length(jrec) > 0 && efrev[jrec[1]] > efrev[jpk] + 1 &&
        (eftrev[jrec[1]] - efrev[jrec[1]]) < (eftrev[jpk] - efrev[jpk]) - 1,
      if (length(jrec))
        sprintf("ejection fraction %.1f%% at the peak burden on day %d against %.1f%% from tissue alone; %.1f%% against %.1f%% by day %d, as the burden halves",
                efrev[jpk], jpk - 1, eftrev[jpk], efrev[jrec[1]], eftrev[jrec[1]], jrec[1] - 1)
      else "the burden never halves in two years")
PI0 <- PI; PI0[["phi_inf"]] <- 0
PI8 <- PI; PI8[["phi_inf"]] <- 0.8
check("a larger phi_inf gives a deeper nadir, so the parameter is live",
      min(qsp2_lvef(ri$v, ri$d, PI8, sat = sat_i)) <
        min(qsp2_lvef(ri$v, ri$d, PI0, sat = sat_i)) - 5,
      sprintf("phi_inf runs from 0 to 0.8; at this burden the nadir moves from %.1f%% to %.1f%%",
              min(qsp2_lvef(ri$v, ri$d, PI0, sat = sat_i)),
              min(qsp2_lvef(ri$v, ri$d, PI8, sat = sat_i))))

# the grid used by the rest of this section
SG8 <- 10^seq(-3, 3, length.out = 13); TG8 <- 2^seq(-1.5, 1.5, by = 0.75)
TRI8 <- data.frame(mu = -3, sigma = 2, sd_tau = 0.6)
trM8 <- qsp2_pop_traj(PE, qsp2_regimen("nivo_q2w"), SG8, TG8)
trC8 <- qsp2_pop_traj(PE, qsp2_regimen("combo_nivo1_ipi3"), SG8, TG8)
check("every grid point of the test population solved in full", all(trM8$ok) && all(trC8$ok),
      sprintf("%d and %d of %d", sum(trM8$ok), sum(trC8$ok), length(trM8$ok)))
check("T2 is never reached before T1 on the same assay",
      all(is.na(trM8$t2c) | trM8$t2c >= trM8$t1ref, na.rm = TRUE),
      "the day of T2 needs troponin above the limit that day")

CO8 <- list(
  q3w = qsp2_pop_cohort("every 3 weeks, troponin I above 34 ng/L", days = seq(21, 365, by = 21),
                        kin = "TnI", url = 34, limit = 34),
  s0  = qsp2_pop_cohort("clinical diagnosis only", kin = "TnI"),
  wk  = qsp2_pop_cohort("weekly troponin T, doubling to 28 ng/L, reported baseline",
                        days = seq(7, 182, by = 7), kin = "TnT", url = 14, limit = 14,
                        rule = "double", abs_min = 28, base = list(med = 10.5, sdlog = 0.62), H = 182))
prep8 <- qsp2_pop_prep(trM8, PP8)
o8 <- qsp2_pop_obs(trM8, prep8, CO8$q3w)
S8 <- o8$S

# (d) internal consistency of the engine
check("the presentation classes of every grid point have probabilities summing to 1",
      max(abs(rowSums(prep8$WM)[trM8$ok] - 1)) < 1e-12,
      sprintf("largest departure %.1e", max(abs(rowSums(prep8$WM)[trM8$ok] - 1))))
check("severe disease happens only within myocarditis, and the bands add up to it",
      all(S8[, "sev"] <= S8[, "t2"] + 1e-12) && all(S8[, "sev0"] <= S8[, "t2"] + 1e-12) &&
        max(abs(S8[, "sevb40"] + S8[, "sevb30"] + S8[, "sevb0"] - S8[, "sev"])) < 1e-12,
      "T3 is within T2 at every grid point")
check("the day-by-day masses add up to the shares they spread over time",
      max(abs(rowSums(o8$M$dx) - S8[, "dx"])) < 1e-12 &&
        max(abs(rowSums(o8$M$firstpos) - S8[, "pos"])) < 1e-12 &&
        max(abs(rowSums(o8$M$dmg) - S8[, "dmg"])) < 1e-12,
      "diagnoses, first positive results and C3 damage")

# (e) measurement error: a patient presenting when troponin just crosses the limit reads
# positive sometimes, not always, which a noise-free model cannot produce
pn <- PE; pn[["k_sym"]] <- 1; pn[["s_ksym"]] <- 1e-6
on <- qsp2_pop_obs(trM8, qsp2_pop_prep(trM8, qsp2_pop_pp(pn)), CO8$s0)
Wn <- qsp2_pop_W(trM8, TRI8)
fr <- sum(Wn[1, ] * on$S[, "presraised"]) / sum(Wn[1, ] * on$S[, "dxc"])
check("with measurement error, troponin at presentation near the limit is not always positive",
      fr > 0.05 && fr < 0.99, sprintf("%.1f%% of clinical diagnoses read positive", 100 * fr))
# the size of the error is an argument; it has to change the answer
on3 <- qsp2_pop_obs(trM8, qsp2_pop_prep(trM8, qsp2_pop_pp(pn)), CO8$s0, sdlog = 0.3)
fr3 <- sum(Wn[1, ] * on3$S[, "presraised"]) / sum(Wn[1, ] * on3$S[, "dxc"])
check("a larger measurement error changes the share read positive at presentation",
      abs(fr3 - fr) > 0.01, sprintf("%.1f%% at the assay's 9%%, %.1f%% at a log-scale SD of 0.3",
                                    100 * fr, 100 * fr3))

# (f) the baseline: the raised share is the lognormal tail, and s_b moves it
bshare <- function(sb) { pb <- PP8; pb$s_b <- sb; BN <- qsp2_pop_basenodes(CO8$q3w, pb); sum(BN$w[BN$raised]) }
exact <- function(sb) pnorm((log(34) - log(PP8$Tn0 * 34 / QSP2_REF_URL)) / sb, lower.tail = FALSE)
check("the share with a raised baseline is exact and rises with the spread of the baseline",
      abs(bshare(1) - exact(1)) < 1e-12 && bshare(1.5) > bshare(1) && bshare(1) > bshare(0.6),
      sprintf("%.4f, %.4f and %.4f at s_b = 0.6, 1 and 1.5", bshare(0.6), bshare(1), bshare(1.5)))

# (g) background rises: fewer at a higher threshold, and the closed form of the
# expected time above the threshold against numerical integration
pb8 <- qsp2_pop_pbg(c(14, 28, 56, 112), CO8$wk, PP8)
# the integrand is zero below the peak size that reaches 28 ng/L and smooth above it, so the
# integral starts there; integrate()'s default accuracy (about 1e-4) is too coarse for this
z0 <- (log(28) - log(PP8$m_bg)) / PP8$s_bg
Dnum <- integrate(function(z) (log(PP8$m_bg / 28) + PP8$s_bg * z) * dnorm(z), z0, Inf,
                  rel.tol = 1e-12)$value / QSP2_TROPONIN$hs_cTnT$k_el
Dcf <- -log(1 - qsp2_pop_pbg(28, CO8$wk, PP8)) / PP8$lam_bg
check("a higher threshold catches fewer background rises, and the closed form is exact",
      all(diff(pb8) < 0) && abs(Dnum - Dcf) < 1e-8 * Dnum,
      sprintf("per-test probability %s at 14, 28, 56 and 112 ng/L; time above 28 ng/L %.6f days in closed form, %.6f by integration",
              paste(signif(pb8, 3), collapse = ", "), Dcf, Dnum))

# (h) the spread of the symptom threshold changes detection before presentation
o1_at <- function(s) {
  pk <- PP8; pk$s_ksym <- s
  o <- qsp2_pop_obs(trM8, qsp2_pop_prep(trM8, pk), CO8$q3w)
  sum(Wn[1, ] * o$S[, "o1num"]) / sum(Wn[1, ] * o$S[, "o1den"])
}
o1s <- c(o1_at(0), o1_at(0.8))
check("the spread of the symptom threshold moves the share detected before presentation",
      abs(o1s[1] - o1s[2]) > 0.01, sprintf("O1 %.3f with no spread, %.3f with a spread of 0.8", o1s[1], o1s[2]))

# (i) the baseline LVEF: severity given a nadir is the integral over the baseline LVEF
e8 <- QSP2_MYO$EF0 + QSP2_EF_SD * qnorm((1:200000 - 0.5) / 200000)
gc8 <- qsp2_g_crit(e8)
sev_cf <- vapply(c(0.6, 0.75, 0.79), function(g) .sev_c(g, 0, 1)[1, "t3"], 0)
sev_bf <- vapply(c(0.6, 0.75, 0.79), function(g) mean(gc8 > g), 0)
check("severity given the nadir equals the integral over baseline LVEF, and is not a step",
      max(abs(sev_cf - sev_bf)) < 1e-4 && all(sev_cf > 0.01 & sev_cf < 0.9999),
      sprintf("P(T3) %s at G = 0.6, 0.75 and 0.79 (brute force %s)",
              paste(sprintf("%.4f", sev_cf), collapse = ", "), paste(sprintf("%.4f", sev_bf), collapse = ", ")))
check("a patient starting at 60% needs exactly a fall below 50%",
      abs(qsp2_g_crit(60) - 0.8) < 1e-12, sprintf("critical G %.3f", qsp2_g_crit(60)))

# (i2) the two troponin statistics are different quantities: C6 reports each patient's highest
# value during therapy, and the engine also keeps the value at the time of diagnosis
oT <- qsp2_pop_obs(trC8, qsp2_pop_prep(trC8, PP8), qsp2_pop_cohort("S0T", kin = "TnT", url = 14))
medh <- function(h) { cs <- cumsum(h); if (sum(h) <= 0) NA_real_ else
  10^QSP2_TROP_BINS[which(cs >= 0.5 * sum(h))[1]] }
mt <- medh(colSums(oT$trop)); mp <- medh(colSums(oT$troppk))
check("the highest troponin of the episode is above the value at presentation",
      is.finite(mt) && is.finite(mp) && mp >= mt,
      sprintf("median %.0f ng/L at presentation against %.0f ng/L at its highest", mt, mp))

# (i3) assay conversion.  Changing a decision threshold at a fixed assay must change detection,
# and the conversion between troponin I assays must do what it states.  Both are tested by effect.
# The hand-set vector above kills so hard (k_kill 3e-2) that every patient, at every node of the
# grid, has a positive result whatever the assay or the background, so it cannot show these
# mechanisms.  A mild version, with k_kill at the bottom of its prior, has almost no drug-related
# injury, so its positives come from baselines and background rises, which is where the assay,
# the background and the exposure act.
PEL <- PE; PEL[["k_kill"]] <- 1e-5
trML <- qsp2_pop_traj(PEL, qsp2_regimen("nivo_q2w"), SG8, TG8)
trCL <- qsp2_pop_traj(PEL, qsp2_regimen("combo_nivo1_ipi3"), SG8, TG8)
prepL <- qsp2_pop_prep(trML, PP8)
obsA <- function(...) qsp2_pop_obs(trML, prepL, qsp2_pop_cohort("assay", days = seq(21, 365, by = 21), ...))
WA <- qsp2_pop_W(trML, TRI8)[1, ]
posA <- function(o) sum(WA * o$S[, "pos"])
aRef <- obsA(kin = "TnI", url = 26.2, limit = 26.2); a55 <- obsA(kin = "TnI", url = 55, limit = 55)
aHi <- obsA(kin = "TnI", url = 26.2, limit = 52.4); aNone <- obsA(kin = "TnI", url = 55, limit = 55, conv = "none")
check("on one troponin I assay, a higher decision threshold finds fewer positives",
      posA(aHi) < posA(aRef) - 1e-4,
      sprintf("share positive %.4f at 26.2 ng/L and %.4f at 52.4 ng/L on the reference assay", posA(aRef), posA(aHi)))
check("another troponin I assay at its own 99th percentile behaves as the reference at its own",
      max(abs(a55$S - aRef$S)) < 1e-10,
      sprintf("largest difference over every statistic and node %.1e; the alternative with no conversion gives %.4f against %.4f",
              max(abs(a55$S - aRef$S)), posA(aNone), posA(aRef)))
check("the alternative with no conversion is live", abs(posA(aNone) - posA(aRef)) > 1e-4,
      sprintf("%.4f against %.4f", posA(aNone), posA(aRef)))

# (j) the engine against the patient-by-patient Monte Carlo.  Besides the basic cohorts, four
# are added: background following the baseline at the top of its prior (g_bg = 2), once
# on a cohort with a reported raised-baseline stratum; a testing exposure shorter than the
# calendar; and a cohort built like C5 (another troponin I assay at its own limit, every
# on-treatment value counted, the testing exposure).  Chronology (no test after a
# clinical diagnosis) is in every cohort with tests.
PP8g <- PP8; PP8g$g_bg <- 2
co_expo <- qsp2_pop_cohort("every 3 weeks, mean 4 tests after baseline", days = seq(21, 365, by = 21),
                           kin = "TnI", url = 34, limit = 34, expo = list(mean = 4, max = 34))
co_c5 <- qsp2_pop_cohort("like C5: 55 ng/L, threshold rule, exposure", days = seq(21, 365, by = 21),
                         kin = "TnI", url = 55, limit = 55, rule = "threshold",
                         expo = list(mean = 7.74, max = 34))
# each added cohort runs on the sick hand-set vector, where diagnosis, severity and chronology are
# exercised, and on the mild one, where the background, the baseline and the exposure are
CO8x <- list(
  q3w_gbg2      = list(co = CO8$q3w, pp = PP8g, tr = trM8, tri = TRI8),
  wk_gbg2       = list(co = CO8$wk,  pp = PP8g, tr = trC8, tri = TRI8),
  q3w_expo      = list(co = co_expo, pp = PP8,  tr = trM8, tri = TRI8),
  c5like        = list(co = co_c5,   pp = PP8g, tr = trC8, tri = TRI8),
  q3w_mild      = list(co = CO8$q3w, pp = PP8,  tr = trML, tri = TRI8),
  q3w_gbg2_mild = list(co = CO8$q3w, pp = PP8g, tr = trML, tri = TRI8),
  wk_gbg2_mild  = list(co = CO8$wk,  pp = PP8g, tr = trCL, tri = TRI8),
  q3w_expo_mild = list(co = co_expo, pp = PP8,  tr = trML, tri = TRI8),
  c5like_mild   = list(co = co_c5,   pp = PP8g, tr = trCL, tri = TRI8))
MC_N <- as.integer(Sys.getenv("QSP2_B4_MC_N", "30000"))
cmp8 <- NULL
runs8 <- c(lapply(setNames(names(CO8), names(CO8)), function(nm)
  list(co = CO8[[nm]], pp = PP8, tr = if (nm == "q3w") trM8 else trC8, tri = TRI8)), CO8x)
for (nm in names(runs8)) {
  trk <- runs8[[nm]]$tr; trik <- runs8[[nm]]$tri
  e <- qsp2_engine_summary(trk, runs8[[nm]]$pp, runs8[[nm]]$co, trik)
  mcs <- qsp2_mc_summary(qsp2_mc_cohort(trk, runs8[[nm]]$pp, runs8[[nm]]$co, trik, MC_N, seed = 8))
  stat <- intersect(names(e), rownames(mcs))
  z <- data.table(cohort = nm, stat = stat, engine = unname(e[stat]), mc = mcs[stat, 1], se = mcs[stat, 2])
  z[, time := grepl("tdx|tfp|tdmg", stat)]
  z[, tol := ifelse(time, pmax(2, 0.1 * abs(mc)), pmax(3 * ifelse(is.na(se), 0, se), 0.1 * abs(mc), 0.002))]
  z[, ok := (is.na(engine) & is.na(mc)) | (is.finite(engine) & is.finite(mc) & abs(engine - mc) <= tol)]
  cmp8 <- rbind(cmp8, z)
}
fwrite(cmp8, file.path(OUT_DIR, "B04_engine_vs_mc.csv"))
print(cmp8[, .(cohort, stat, engine = signif(engine, 4), mc = signif(mc, 4), tol = signif(tol, 3), ok)])
check("the population engine agrees with a Monte Carlo of individual patients",
      all(cmp8$ok), sprintf("%d of %d statistics within tolerance in %d cohorts; worst: %s",
                            sum(cmp8$ok), nrow(cmp8), uniqueN(cmp8$cohort),
                            cmp8[which.max(abs(engine - mc) / tol), paste(cohort, stat)]))
# the added cohorts must exercise what they were added for, or agreeing on them proves nothing
mcpos <- function(nm) cmp8[cohort == nm & stat == "pos", mc]
check("the comparison exercises the baseline-dependent background",
      abs(mcpos("q3w_gbg2_mild") - mcpos("q3w_mild")) > 3 * cmp8[cohort == "q3w_mild" & stat == "pos", se] + 0.005 &&
        mcpos("q3w_mild") < 0.95,
      sprintf("share with a positive result %.3f at g_bg = 0 and %.3f at g_bg = 2 (Monte Carlo, mild vector)",
              mcpos("q3w_mild"), mcpos("q3w_gbg2_mild")))
check("the comparison exercises the testing exposure",
      mcpos("q3w_expo_mild") < mcpos("q3w_mild") - 0.005 && mcpos("q3w_mild") < 0.95,
      sprintf("share with a positive result %.3f with every test, %.3f with a mean of 4 tests (Monte Carlo, mild vector)",
              mcpos("q3w_mild"), mcpos("q3w_expo_mild")))

# (k) the C6 regimen-mix helper reads the statistic the calibration fits: at the cohort's own
# 10% mix it must return exactly the C6 statistic of qsp2_cal_stats.
source(file.path(QSP2_DIR, "R", "lib_lhs.R"), local = TRUE)
drawK <- as.list(PE[intersect(names(QSP2_PRIOR), names(PE))])
runK <- qsp2_cal_run(qsp2_cal_pars(drawK), qsp2_cal_cohorts()["S0T"], SG8)
stK <- qsp2_cal_stats(runK, TRI8)
tmK <- qsp2_cal_tropmix(list(draw = drawK, tri = TRI8), list(sgrid = SG8, mix = c(0, 0.10, 0.3)))
check("the C6 regimen-mix helper returns the C6 statistic the calibration fits",
      is.finite(stK[1, "C6_trop_dx"]) && abs(tmK[2] - stK[1, "C6_trop_dx"]) < 1e-9 * stK[1, "C6_trop_dx"],
      sprintf("%.1f ng/L from the helper at 10%% combination, %.1f from the calibration statistic",
              tmK[2], stK[1, "C6_trop_dx"]))
check("the comparison is not trivial: the test population has severe cases screening prevents",
      cmp8[cohort == "q3w" & stat == "sev0", mc] > cmp8[cohort == "q3w" & stat == "sev", mc] &&
        cmp8[cohort == "q3w" & stat == "O3", mc] < 1 && cmp8[cohort == "q3w" & stat == "O3", mc] > 0,
      sprintf("T3 %.3f without screening, %.3f with it; O3 %.2f",
              cmp8[cohort == "q3w" & stat == "sev0", mc], cmp8[cohort == "q3w" & stat == "sev", mc],
              cmp8[cohort == "q3w" & stat == "O3", mc]))

# =====================================================================================
cat("
9. Regression tests
")
# The same invariants, pinned as testthat tests so that every later block can re-run
# them cheaply before it does anything else.
source(file.path(QSP2_DIR, "R", "run_tests.R"), local = TRUE)
tr <- qsp2_run_tests()
check("the regression test suite passes", isTRUE(tr$ok),
      sprintf("%d tests, %d failed, %d errors", tr$n, tr$failed, tr$errors))

# =====================================================================================
res <- copy(.checks)
fwrite(res, file.path(OUT_DIR, "B04_verification.csv"))
cat(sprintf("\n%d checks recorded in output/B04_verification.csv\n", nrow(res)))
save_json(list(block = "B04",
               n_states = length(y0), n_parameters = length(pars),
               rhs_worst_relative_difference = worst, random_states = nrand,
               half_life_nivolumab_days = unname(h_n[["beta"]]),
               half_life_ipilimumab_days = unname(h_i[["beta"]]),
               half_life_pembrolizumab_days = unname(h_p[["beta"]]),
               spec_sha256 = sha256_file(file.path(DOC_DIR, "model_specification.md")),
               tests = tr$n, tests_failed = tr$failed,
               engine_vs_mc_within_tolerance = sum(cmp8$ok), engine_vs_mc_compared = nrow(cmp8),
               tnt_qss_worst_error = max(qss[, "err"])),
          "B04_implementation.json")
