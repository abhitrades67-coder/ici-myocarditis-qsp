# ---------------------------------------------------------------------------
# lib_cal.R  -  the calibration cohorts, targets and statistics, shared by blocks B6 to B8.
#
# This file turns each calibration statistic into the population quantity the engine of
# lib_pop.R computes, so every block reads the same definitions.
#
# Definitions only: sourcing this file runs nothing.
# ---------------------------------------------------------------------------

QSP2_SGRID_CAL  <- 10^seq(-4, 4, length.out = 25)     # calibration search
QSP2_SGRID_FULL <- 10^seq(-4, 4, length.out = 41)     # finer grid, for convergence checks

# The population search grid: the lognormal of S (mu, sigma) and the spread of tempo, whose
# prior is uniform from 0.05 to 1.5.  sigma runs to 4.4, beyond which the central 95% of S no
# longer fits on the S grid (ln 1e8 = 18.4 = 2 x 1.96 x 4.7) and more spread only moves mass
# into the end cells; mu runs down to -12, so that most of a population may sit below the S
# grid, which the end-cell rule (qsp2_cal_edges) allows only where the lowest node is free of
# drug-related injury.
QSP2_TRI_GRID <- expand.grid(mu = seq(-12, 4, by = 0.5), sigma = seq(0.4, 4.4, by = 0.25),
                             sd_tau = c(0.05, 0.15, 0.3, 0.5, 0.75, 1, 1.25, 1.5))

# Trajectories are solved to 730 days, because Cheng (C5) follows patients for 24 months and
# van den Berg (C6) for a mean of 1.6 years; cohorts S0 and C3 keep their own shorter follow-up.
QSP2_CAL_H <- 730

# Testing exposure: the number of tests after the baseline, geometric on 1 to 34 with mean
# 7.74.  van den Berg 2024: "The frequency of hs-TnT measurements per individual patient
# spanned a range from 2 to 35, resulting in a mean value of 8.74 measurements per patient", with
# testing "at each treatment visit where new dose was administered", so 8.74 includes the
# baseline.  It is the only reported treatment exposure among the screening cohorts, and is
# applied to C4 and C5 as well, as an assumption.
QSP2_EXPO <- list(mean = 7.74, max = 34)

# The C1 case file: FAERS myocarditis cases with onset counted from the earliest start of a
# suspect checkpoint inhibitor, with the precision of each date, built from the FAERS quarterly
# files by R/faers_c1_cases.R.  C1 uses only cases with both dates recorded to the day.
QSP2_C1_FILE <- function() file.path(DATA_DIR, "c1_cases_precision.csv")

# The prior of the fitted parameters.  Each entry is c(low, high, kind).
QSP2_PRIOR <- list(
  lam_p = c(1e2, 1e6, "log"), gC = c(0.1, 100, "log"), gPp = c(0.1, 100, "log"),
  rho = c(0.05, 1.5, "unif"), gP = c(0.1, 100, "log"), gF = c(0.1, 10, "log"),
  k_tr = c(0.05, 1, "unif"), NT = c(2, 8, "int"), phi = c(0.01, 0.6, "log"),
  pi_h = c(0, 0.15, "unif"), del_h = c(0.05, 0.5, "unif"),
  k_kill = c(1e-5, 1e-2, "log"), Kh = c(1e6, 1e9, "log"),
  k_lic = c(0.02, 2, "log"), K_lic = c(1e4, 1e7, "log"),
  wL = c(0, 5, "unif"), wM = c(0, 5, "unif"), eps_x = c(0.001, 0.5, "log"),
  k_x = c(1, 200, "log"), r_rep = c(0.001, 0.5, "log"), r_scar = c(0.001, 0.5, "log"),
  fc = c(0.05, 0.5, "unif"), k_rel = c(0.02, 0.5, "log"),
  eta = c(0, 1, "unif"), pexp = c(0.5, 4, "unif"),
  EC50pd1 = c(0.05, 5, "log"), EC50ctla4 = c(0.3, 6, "log"),
  f_dx = c(0.03, 1, "log"), k_sym = c(1.5, 50, "log"),
  s_ksym = c(0, 1, "unif"), L_tx = c(1, 14, "unif"),
  s_b = c(0.3, 2, "unif"), lam_bg = c(1e-3, 1, "log"), m_bg = c(10, 200, "log"),
  # the greatest fraction of function inflammation can take away
  phi_inf = c(0, 0.8, "unif"),
  # zero switches the term off, so the prior widens the model rather than moving it
  g_bg = c(0, 2, "unif"))

#' The calibration cohorts, each with its own calendar, assay, threshold, rule, baseline,
#' follow-up and regimen mix.
qsp2_cal_cohorts <- function() list(
  # C1 and C2: diagnosis at clinical presentation only, pure regimens
  S0  = qsp2_pop_cohort("S0, hs-cTnI (C1, C2)", kin = "TnI", url = QSP2_REF_URL, H = 365),
  # C6, Pereyra Pietri: troponin T at the diagnosis of clinically recognised myocarditis.  The
  # source does not report its regimen mix; 10% combination is assumed and B8 reports how
  # far the statistic moves with it.
  S0T = qsp2_pop_cohort("S0, hs-cTnT (C6 Pereyra Pietri)", kin = "TnT", url = QSP2_TNT_URL,
                        H = 365, w_combo = 0.10),
  # C3, Nishikawa: monthly for 6 months, Abbott hs-cTnI above 26.2 ng/L, 5% combination
  C3 = qsp2_pop_cohort("C3 Nishikawa", days = seq(30, 180, by = 30), kin = "TnI",
                       url = QSP2_REF_URL, limit = QSP2_REF_URL, H = 180, w_combo = 0.05),
  # C4, Toublanc: before each infusion, "every 2 to 4 weeks ... led by anti-PD(L)1 administration
  # schedule", hs-cTnI above 34 ng/L (the source's own 99th percentile).  The calendar follows
  # the drug mix of Table 1, pembrolizumab 198, atezolizumab 17 and nivolumab
  # plus ipilimumab 7 every 3 weeks (222 of 298), nivolumab 34 and durvalumab 42 every 2 or 4
  # weeks, split equally because the source does not say which; nivolumab plus ipilimumab in 7.
  # Follow-up to the end of treatment, so the exposure model sets how long testing lasts.
  C4 = qsp2_pop_cohort_mix("C4 Toublanc", w = c(38, 222, 38) / 298, parts = lapply(c(14, 21, 28),
         function(iv) qsp2_pop_cohort(sprintf("C4 Toublanc, every %d days", iv),
                                      days = seq(iv, QSP2_CAL_H, by = iv), kin = "TnI", url = 34,
                                      limit = 34, H = QSP2_CAL_H, w_combo = 7 / 298, expo = QSP2_EXPO))),
  # C5, Cheng: before each dose, Siemens cTnI at or above 55 ng/L (the assay's own 99th
  # percentile), anti-CTLA-4 plus PD-1 in 56 of 428, followed for 24 months.
  # "The initial elevated cTnI was defined as the first elevated value during ICI therapy", and
  # 13 of the 42 counted had a known raised baseline, so every on-treatment value at or above the
  # limit counts whatever the baseline (rule "threshold"), and the baseline sample itself does
  # not.  The source gives no schedule by drug, so every 3 weeks is assumed.
  C5 = qsp2_pop_cohort("C5 Cheng", days = seq(21, QSP2_CAL_H, by = 21), kin = "TnI", url = 55,
                       limit = 55, rule = "threshold", H = QSP2_CAL_H, w_combo = 56 / 428,
                       expo = QSP2_EXPO),
  # C6, van den Berg: hs-TnT at each treatment visit "every 3-4 weeks", 2 to 35 measurements per
  # patient (mean 8.74), mean follow-up 1.60 years; a significant rise is a doubling to at least
  # 28 ng/L; baseline median 10.5 ng/L (IQR 6.5 to 15); PD-(L)1 plus CTLA-4 in 58 of 164.  Half
  # the patients are tested every 3 weeks and half every 4, each with the reported exposure.
  C6 = qsp2_pop_cohort_mix("C6 van den Berg", w = c(0.5, 0.5), parts = lapply(c(21, 28),
         function(iv) qsp2_pop_cohort(sprintf("C6 van den Berg, every %d days", iv),
                                      days = seq(iv, QSP2_CAL_H, by = iv), kin = "TnT",
                                      url = QSP2_TNT_URL, limit = QSP2_TNT_URL, rule = "double",
                                      abs_min = 28,
                                      base = list(med = 10.5, sdlog = log(15 / 6.5) / (2 * qnorm(0.75))),
                                      H = QSP2_CAL_H, w_combo = 58 / 164, expo = QSP2_EXPO))))

#' Bootstrap variances of the C1 statistics from the case file.
# C1 is defined over a 12-month horizon (data/targets_calibration.csv), which is also how long
# cohort S0 is simulated, so a reported onset after day 365 is left out: the model cannot place
# a case there.  Reports from Japan are left out.  Only cases whose start and event dates are
# both recorded to the day are used: a date given only as a month or a year is completed to its
# first day in FAERS processing, which creates spurious zero-day onsets.  A genuine same-day case
# with both dates complete stays in.  `complete = FALSE` keeps the completed dates as well.
qsp2_c1_bootstrap <- function(file = QSP2_C1_FILE(), B = 20000, seed = 20260928, horizon = 365,
                              complete = TRUE) {
  x <- read.csv(file, stringsAsFactors = FALSE)
  x$occr_country[is.na(x$occr_country)] <- ""      # reports with no country are kept
  x <- x[x$occr_country != "JP", ]
  n_all <- c(mono = sum(x$regimen == "PD-1 monotherapy" & x$onset_days <= horizon),
             combo = sum(x$regimen == "COMBO (PD-[L]1 + CTLA-4)" & x$onset_days <= horizon))
  if (complete) x <- x[x$start_prec == "day" & x$event_prec == "day", ]
  mo_all <- x$onset_days[x$regimen == "PD-1 monotherapy"]
  co_all <- x$onset_days[x$regimen == "COMBO (PD-[L]1 + CTLA-4)"]
  mo <- mo_all[mo_all <= horizon]
  co <- co_all[co_all <= horizon]
  set.seed(seed)
  bm <- replicate(B, quantile(sample(mo, replace = TRUE), c(.25, .5, .75), names = FALSE))
  bc <- replicate(B, quantile(sample(co, replace = TRUE), c(.25, .5, .75), names = FALSE))
  list(n_mono = length(mo), n_combo = length(co), horizon = horizon, complete = complete,
       n_mono_any_date = unname(n_all["mono"]), n_combo_any_date = unname(n_all["combo"]),
       beyond_mono = sum(mo_all > horizon), beyond_combo = sum(co_all > horizon),
       q_mono = quantile(mo, c(.25, .5, .75), names = FALSE),
       q_combo = quantile(co, c(.25, .5, .75), names = FALSE),
       sd_mono = apply(bm, 1, sd), sd_combo = apply(bc, 1, sd),
       sd_ratio = sd(bm[2, ] / bc[2, ]))
}

#' The calibration targets: observed value, observational variance, and whether the value
#' is a point or a lower bound.  The model discrepancy is 10% of the observed value for every
#' target.  `r2` marks the incidence values requirement R2 uses.
qsp2_cal_targets <- function(c1) {
  bin <- function(k, n) (k / n) * (1 - k / n) / n
  q25 <- function(v) (0.25 * v)^2      # a median or quartile reported without an interval
  # Pereyra Pietri: median 364 ng/L, IQR 1096, in 59 patients.  Read as lognormal, the IQR
  # gives a log-scale SD s = asinh(1096 / (2 * 364)) / qnorm(0.75), and the sampling SD of a
  # median is m * s * sqrt(2 * pi) / (2 * sqrt(n)).
  s_pp <- asinh(1096 / (2 * 364)) / qnorm(0.75)
  sd_pp <- 364 * s_pp * sqrt(2 * pi) / (2 * sqrt(59))
  t <- rbind(
    # read from the case file: within the 12-month horizon, cases with both dates complete to
    # the day.  The values are taken from c1, not written here by hand, so the definition and
    # its numbers cannot drift apart.
    data.frame(id = "C1_q25_combo", obs = c1$q_combo[1], v_obs = c1$sd_combo[1]^2, kind = "point", r2 = FALSE),
    data.frame(id = "C1_q50_combo", obs = c1$q_combo[2], v_obs = c1$sd_combo[2]^2, kind = "point", r2 = FALSE),
    data.frame(id = "C1_q75_combo", obs = c1$q_combo[3], v_obs = c1$sd_combo[3]^2, kind = "point", r2 = FALSE),
    data.frame(id = "C1_q25_mono",  obs = c1$q_mono[1],  v_obs = c1$sd_mono[1]^2,  kind = "point", r2 = FALSE),
    data.frame(id = "C1_q50_mono",  obs = c1$q_mono[2],  v_obs = c1$sd_mono[2]^2,  kind = "point", r2 = FALSE),
    data.frame(id = "C1_q75_mono",  obs = c1$q_mono[3],  v_obs = c1$sd_mono[3]^2,  kind = "point", r2 = FALSE),
    data.frame(id = "C1_ratio",     obs = c1$q_mono[2] / c1$q_combo[2], v_obs = c1$sd_ratio^2,
               kind = "point", r2 = FALSE),
    data.frame(id = "C2_p_combo",   obs = 0.0027, v_obs = bin(0.0027 * 2686, 2686), kind = "lower", r2 = TRUE),
    data.frame(id = "C2_p_mono",    obs = 0.0006, v_obs = bin(0.0006 * 17908, 17908), kind = "lower", r2 = TRUE),
    data.frame(id = "C2_ratio",     obs = 4.5,   v_obs = 1.5^2,            kind = "point", r2 = TRUE),
    data.frame(id = "C3_p_damage",  obs = 0.10,  v_obs = bin(10, 100),     kind = "point", r2 = TRUE),
    data.frame(id = "C3_t_median",  obs = 82,    v_obs = q25(82),          kind = "point", r2 = FALSE),
    data.frame(id = "C4_p_myo",     obs = 5 / 298, v_obs = bin(5, 298),    kind = "point", r2 = TRUE),
    data.frame(id = "C4_t_mean",    obs = 117,   v_obs = 51^2 / 5,         kind = "point", r2 = FALSE),
    data.frame(id = "C4_p_screen",  obs = 0.478, v_obs = 0,                kind = "lower", r2 = FALSE),
    data.frame(id = "C4_ppv",       obs = 0.20,  v_obs = bin(5, 25),       kind = "point", r2 = FALSE),
    data.frame(id = "C4_p_baseline", obs = 15 / 298, v_obs = bin(15, 298), kind = "point", r2 = FALSE),
    data.frame(id = "C5_p_pos",     obs = 42 / 428, v_obs = bin(42, 428),  kind = "point", r2 = TRUE),
    data.frame(id = "C5_q25",       obs = 21.3,  v_obs = q25(21.3),        kind = "point", r2 = FALSE),
    data.frame(id = "C5_q50",       obs = 56,    v_obs = q25(56),          kind = "point", r2 = FALSE),
    data.frame(id = "C5_q75",       obs = 96.3,  v_obs = q25(96.3),        kind = "point", r2 = FALSE),
    data.frame(id = "C5_ppv",       obs = 6 / 42, v_obs = bin(6, 42),      kind = "point", r2 = TRUE),
    data.frame(id = "C5_t_myo",     obs = 82,    v_obs = q25(82),          kind = "point", r2 = FALSE),
    data.frame(id = "C6_p_rise_nodx", obs = 18 / 164, v_obs = bin(18, 164), kind = "point", r2 = FALSE),
    data.frame(id = "C6_trop_dx",   obs = 364,   v_obs = sd_pp^2,          kind = "point", r2 = FALSE),
    data.frame(id = "C6_m_bg",      obs = 41.9,  v_obs = q25(41.9),        kind = "point", r2 = FALSE))
  t$v_md <- (0.10 * t$obs)^2
  t
}

#' Targets kept as sensitivity analyses only, with the same fields.
#' C5_p_presented read Cheng's "83.3% were classified as grade 3 or higher" as the share of
#' diagnosed cases that had presented clinically before screening found them.  Grade is not that
#' observable: the source reports that 4 of the 5 symptomatic cases "had not initially sought
#' medical attention for their symptoms; rather, they were first flagged by elevated cTnI".
#' It is not fitted; block B8 records whether each population would also meet it.
qsp2_cal_targets_sens <- function() {
  t <- data.frame(id = "C5_p_presented", obs = 5 / 6, v_obs = (5 / 6) * (1 / 6) / 6,
                  kind = "point", r2 = FALSE)
  t$v_md <- (0.10 * t$obs)^2
  t
}

#' A full parameter vector from a draw of the prior (a named list or one-row table).
qsp2_cal_pars <- function(draw) {
  p <- qsp2_default_pars(assay = "hs_cTnI")
  for (nm in intersect(names(draw), names(p))) p[[nm]] <- as.numeric(draw[[nm]])
  p[["NT"]] <- round(p[["NT"]])
  p
}

#' Solve both regimens on the grid and observe every cohort.  `regimens` names the monotherapy
#' and combination calendars.
# The trajectories are solved as far as the longest follow-up among the cohorts (730 days, with
# the regimen's doses continuing that long); a cohort is observed over its own follow-up.
# Monotherapy is solved as nivolumab every 2 weeks for every cohort; the calendars of the tests
# follow each cohort's drug mix.  Pembrolizumab every 3 weeks gives the same T2 incidence and
# median day of T2 as nivolumab every 2 weeks, because PD-1 occupancy is saturated under both.
qsp2_cal_run <- function(pars, cohorts = qsp2_cal_cohorts(), sgrid = QSP2_SGRID_CAL,
                         regimens = c(mono = "nivo_q2w", combo = "combo_nivo1_ipi3"),
                         no_treatment = FALSE) {
  HT <- max(vapply(cohorts, function(co) co$H, numeric(1)))
  tr <- lapply(regimens, function(r) qsp2_pop_traj(pars, qsp2_regimen(r, horizon = HT), sgrid, H = HT))
  pp <- qsp2_pop_pp(pars)
  Hs <- sort(unique(vapply(cohorts, function(co) min(co$H, HT), numeric(1))))
  prep <- lapply(tr, function(t) setNames(lapply(Hs, function(h) qsp2_pop_prep(t, pp, H = h)),
                                          as.character(Hs)))
  obs <- lapply(cohorts, function(co) {
    h <- as.character(min(co$H, HT))
    setNames(lapply(names(tr), function(r)
      qsp2_pop_obs_any(tr[[r]], prep[[r]][[h]], co, no_treatment = no_treatment)), names(tr))
  })
  list(pars = pars, tr = tr, obs = obs, cohorts = cohorts)
}

#' The end-cell rule, for every population of `tri`.  The cell
#' masses of qsp2_pop_cellmass run to zero and to infinity at the ends of the S grid, so a
#' population wider than the grid puts its tails on the two end nodes.  At the low end that is
#' harmless if the lowest node has no drug-related injury at all (no T1 at any tempo in either
#' regimen): patients with still less precursor supply have none either, and the lowest node
#' stands for them correctly.  At the high end there is no such guarantee, so the highest cell
#' may hold at most 1% of the population.  The tempo grid keeps its end cells under 1%.
qsp2_cal_edges <- function(run, tri) {
  low <- vapply(run$tr, function(t) all(is.na(t$t1ref[t$nodes$i_s == 1L])), logical(1))
  sg <- run$tr[[1]]$sgrid
  m <- vapply(seq_len(nrow(tri)), function(i) {
    w <- qsp2_pop_cellmass(sg, tri$mu[i], tri$sigma[i]); c(w[1], w[length(w)]) }, numeric(2))
  cbind(edge_lo = m[1, ], edge_hi = m[2, ], edge_low_free = all(low),
        edge_ok = as.numeric((m[1, ] <= 0.01 | all(low)) & m[2, ] <= 0.01))
}

#' Every calibration statistic, and the quantities requirements R1, R3 and R4 read, for
#' every population on the search grid: a (populations x statistics) matrix.
qsp2_cal_stats <- function(run, tri) {
  WM <- qsp2_pop_W(run$tr$mono, tri); WC <- qsp2_pop_W(run$tr$combo, tri)
  O <- run$obs; CO <- run$cohorts
  n <- nrow(tri)
  wof <- function(c, w) if (is.null(w)) CO[[c]]$w_combo else w
  # a block may run only some cohorts (B6 needs S0 and C3 to C5); the rest come back NA
  sh <- function(c, what, w = NULL) {
    if (!c %in% names(O)) return(rep(NA_real_, n))
    w <- wof(c, w)
    a <- if (w < 1) as.vector(WM %*% O[[c]]$mono$S[, what]) else 0
    b <- if (w > 0) as.vector(WC %*% O[[c]]$combo$S[, what]) else 0
    (1 - w) * a + w * b
  }
  cm <- function(c, mass, w = NULL) {
    if (!c %in% names(O)) return(matrix(0, n, 1))
    w <- wof(c, w)
    a <- if (w < 1) WM %*% matrixStats::rowCumsums(O[[c]]$mono$M[[mass]]) else 0
    b <- if (w > 0) WC %*% matrixStats::rowCumsums(O[[c]]$combo$M[[mass]]) else 0
    (1 - w) * a + w * b
  }
  tq <- function(c, mass, probs, w = NULL) {
    if (!c %in% names(O)) return(matrix(NA_real_, n, length(probs)))
    qsp2_pop_tq_cum(cm(c, mass, w), probs)
  }
  rat <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)
  q1c <- tq("S0", "dxc", c(.25, .5, .75), 1); q1m <- tq("S0", "dxc", c(.25, .5, .75), 0)
  p2c <- sh("S0", "dxc", 1); p2m <- sh("S0", "dxc", 0)
  q5 <- tq("C5", "firstpos", c(.25, .5, .75))
  # an alternative reading of C1, recorded but not used for acceptance: reported onset as a
  # mixture of clinical diagnosis (three quarters) and diagnosis by screening before each dose
  # (a quarter)
  mixq <- function(w) {
    if (!all(c("S0", "C5") %in% names(O))) return(matrix(NA_real_, n, 3))
    a <- cm("S0", "dxc", w); b <- cm("C5", "dx", w)
    a <- a / pmax(a[, ncol(a)], 1e-300); b <- b / pmax(b[, ncol(b)], 1e-300)
    qsp2_pop_tq_cum(0.75 * a + 0.25 * b, c(.25, .5, .75))
  }
  x1c <- mixq(1); x1m <- mixq(0)
  # C6's troponin statistic is each patient's highest value during therapy, not the value at the
  # moment of presentation: the source's Methods say "Highest hs-TnT value during ICI therapy was
  # considered".  The value at presentation is a different quantity and is biased low.
  trop6 <- if ("S0T" %in% names(O)) {
    wT <- CO$S0T$w_combo
    qsp2_pop_tropmed_hist((1 - wT) * (WM %*% O$S0T$mono$troppk) + wT * (WC %*% O$S0T$combo$troppk))
  } else rep(NA_real_, n)
  # R1: among T2 cases, the first day of T1 on the reference assay
  r1 <- function(tr, W) {
    t2n <- as.numeric(!is.na(tr$t2c) & tr$t2c <= 365)
    den <- as.vector(W %*% t2n)
    p7 <- rat(as.vector(W %*% (t2n * (!is.na(tr$t1ref) & tr$t1ref >= 7))), den)
    cbind(p7 = p7, med = qsp2_pop_nodeq(W, tr$t1ref, t2n, 0.5)[, 1], t2 = den)
  }
  R1m <- r1(run$tr$mono, WM); R1c <- r1(run$tr$combo, WC)
  st <- cbind(
    C1_q25_combo = q1c[, 1], C1_q50_combo = q1c[, 2], C1_q75_combo = q1c[, 3],
    C1_q25_mono = q1m[, 1], C1_q50_mono = q1m[, 2], C1_q75_mono = q1m[, 3],
    C1_ratio = rat(q1m[, 2], q1c[, 2]),
    C2_p_combo = p2c, C2_p_mono = p2m, C2_ratio = rat(p2c, p2m),
    C1mix_q25_combo = x1c[, 1], C1mix_q50_combo = x1c[, 2], C1mix_q75_combo = x1c[, 3],
    C1mix_q25_mono = x1m[, 1], C1mix_q50_mono = x1m[, 2], C1mix_q75_mono = x1m[, 3],
    C1mix_ratio = rat(x1m[, 2], x1c[, 2]),
    # an alternative reading of C2, recorded but not used for acceptance: C2 read as T3 under S0
    C2T3_p_combo = sh("S0", "sev", 1), C2T3_p_mono = sh("S0", "sev", 0),
    C2T3_ratio = rat(sh("S0", "sev", 1), sh("S0", "sev", 0)),
    C3_p_damage = sh("C3", "dmg"), C3_t_median = tq("C3", "dmg", 0.5)[, 1],
    C4_p_myo = sh("C4", "dx"), C4_t_mean = rat(sh("C4", "m1dx"), sh("C4", "dx")),
    C4_p_screen = rat(sh("C4", "dxs"), sh("C4", "dx")),
    C4_ppv = rat(sh("C4", "dxpos"), sh("C4", "pos")), C4_p_baseline = sh("C4", "basep"),
    C5_p_pos = sh("C5", "pos"), C5_q25 = q5[, 1], C5_q50 = q5[, 2], C5_q75 = q5[, 3],
    C5_ppv = rat(sh("C5", "dxpos"), sh("C5", "pos")),
    C5_t_myo = tq("C5", "firstposdx", 0.5)[, 1],
    C5_p_presented = rat(sh("C5", "dxpres"), sh("C5", "dx")),
    C6_p_rise_nodx = sh("C6", "pos") - sh("C6", "dxpos"), C6_trop_dx = trop6,
    C6_m_bg = rep(run$pars[["m_bg"]], nrow(tri)),
    # requirements, on monotherapy (the cohorts that set R1's range are 87% to 98%
    # monotherapy) and, for the record, on combination
    R1_p7_mono = R1m[, "p7"], R1_med_mono = R1m[, "med"],
    R1_p7_combo = R1c[, "p7"], R1_med_combo = R1c[, "med"],
    R3_t3_share_mono = rat(sh("S0", "sev", 0), sh("S0", "t2", 0)),
    R3_b40_mono = rat(sh("S0", "sevb40", 0), sh("S0", "sev", 0)),
    R3_b30_mono = rat(sh("S0", "sevb30", 0), sh("S0", "sev", 0)),
    R3_b0_mono = rat(sh("S0", "sevb0", 0), sh("S0", "sev", 0)),
    R3_t3_share_combo = rat(sh("S0", "sev", 1), sh("S0", "t2", 1)),
    R4_t3_mono = tq("S0", "dxcsev", 0.5, 0)[, 1], R4_non_mono = tq("S0", "dxcnon", 0.5, 0)[, 1],
    R4_t3_combo = tq("S0", "dxcsev", 0.5, 1)[, 1], R4_non_combo = tq("S0", "dxcnon", 0.5, 1)[, 1],
    t2_mono = R1m[, "t2"], t2_combo = R1c[, "t2"],
    lost = pmax(attr(WM, "lost"), attr(WC, "lost")))
  st
}

#' Implausibility of each statistic against its target: a (populations x targets) matrix.  A statistic the model cannot produce counts as infinitely implausible.
qsp2_cal_I <- function(st, tg) {
  I <- matrix(NA_real_, nrow(st), nrow(tg), dimnames = list(NULL, tg$id))
  for (j in seq_len(nrow(tg))) {
    f <- st[, tg$id[j]]
    d <- f - tg$obs[j]
    if (tg$kind[j] == "lower") d <- pmin(d, 0)
    I[, j] <- abs(d) / sqrt(tg$v_obs[j] + tg$v_md[j])
    I[!is.finite(f), j] <- Inf
  }
  I
}

#' The largest and second-largest implausibility in each row.
qsp2_top2 <- function(I) {
  m1 <- matrixStats::rowMaxs(I)
  J <- I; J[cbind(seq_len(nrow(I)), max.col(I, ties.method = "first"))] <- -Inf
  cbind(I1 = m1, I2 = if (ncol(I) > 1) matrixStats::rowMaxs(J) else rep(-Inf, nrow(I)))
}

#' Requirements R1 and R3, row by row; both are imposed on accepted sets.
# R1 is the latent phase: the first troponin rise on or after day 7 in at least 95% of T2 cases.
# The cohorts' times to detection are fitted directly through C3_t_median, C4_t_mean, the C5
# quartiles and C5_t_myo.  R1_med_mono, the median day of the first rise, is computed and
# reported.
qsp2_req_R1 <- function(st) with(as.data.frame(st),
  is.finite(R1_p7_mono) & R1_p7_mono >= 0.95)
# R3: fewer than half of T2 cases on monotherapy reach T3 (C4, where all screen-detected
# myocarditis was asymptomatic).  The nadir-LVEF bands are computed and reported by block B6,
# but they are an output and not a requirement.
qsp2_req_R3 <- function(st) with(as.data.frame(st),
  is.finite(R3_t3_share_mono) & R3_t3_share_mono < 0.5)
qsp2_req_R4 <- function(st) with(as.data.frame(st),
  (is.finite(R4_t3_mono) & is.finite(R4_non_mono) & R4_t3_mono < R4_non_mono) |
    (is.finite(R4_t3_combo) & is.finite(R4_non_combo) & R4_t3_combo < R4_non_combo))

#' The worker: one parameter set, every population on the grid.
#'
#' ctx$mode "B6": the population chosen for R2 (the lowest largest implausibility over the
#' R2 values) and R1, R3, R4 there.  ctx$mode "B8": every population that is not implausible
#' over all targets and meets R1 and R3, and the best population for the record.
qsp2_cal_one <- function(job, ctx) {
  res <- tryCatch({
    pars <- qsp2_cal_pars(job$draw)
    if (pars[["rho"]] >= pars[["del_E"]] + pars[["k_tr"]])
      return(list(id = job$id, status = "constraint"))
    t0 <- proc.time()[3]
    run <- qsp2_cal_run(pars, ctx$cohorts, ctx$sgrid)
    st <- cbind(qsp2_cal_stats(run, ctx$tri), qsp2_cal_edges(run, ctx$tri))
    I <- qsp2_cal_I(st, ctx$tg)
    # a population that loses solved mass, or breaks the end-cell rule, is not used
    lost_ok <- st[, "lost"] <= 1e-3 & st[, "edge_ok"] == 1
    out <- list(id = job$id, status = "ok", secs = unname(proc.time()[3] - t0),
                failed_nodes = sum(!run$tr$mono$ok) + sum(!run$tr$combo$ok))
    if (ctx$mode == "B6") {
      I2 <- I[, ctx$tg$id[ctx$tg$r2], drop = FALSE]
      t2i <- qsp2_top2(I2)
      sc <- ifelse(lost_ok, t2i[, 1] + 1e-3 * rowSums(pmin(I2, 1e6)), Inf)
      b <- which.min(sc)
      out$best <- c(tri_row = b, unlist(ctx$tri[b, ]), I2max = unname(t2i[b, 1]), I2second = unname(t2i[b, 2]),
                    st[b, ], setNames(I[b, ], paste0("I_", colnames(I))))
      if (isTRUE(ctx$ipi)) {
        tri_b <- ctx$tri[b, , drop = FALSE]
        tI <- qsp2_pop_traj(pars, qsp2_regimen("ipi_alone"), ctx$sgrid)
        WI <- qsp2_pop_W(tI, tri_b)
        out$best <- c(out$best, R5_t2_ipi = sum(WI[1, ] * (!is.na(tI$t2c) & tI$t2c <= 365)))
      }
    } else {
      t2a <- qsp2_top2(I)
      # the acceptance rule, and for the record the same rule with each source left out in turn,
      # with the C5 grade target added, with C2 read as T3 and with C1 read as a mixture, so
      # that these alternatives can be examined without running the calibration again
      r13 <- lost_ok & qsp2_req_R1(st) & qsp2_req_R3(st)
      passes <- function(Ix) { tt <- qsp2_top2(Ix); tt[, 1] <= 3 & tt[, 2] <= 2.5 }
      ids <- ctx$tg$id; srcs <- sub("_.*", "", ids)
      variants <- c(list(all = ids),
                    setNames(lapply(unique(srcs), function(z) ids[srcs != z]),
                             paste0("no_", unique(srcs))))
      flags <- vapply(variants, function(v) r13 & passes(I[, v, drop = FALSE]), logical(nrow(I)))
      flags <- matrix(flags, nrow(I), dimnames = list(NULL, names(variants)))
      # the C5 grade target is a sensitivity analysis; whether a population would also pass
      # with it added is recorded
      Is <- if (!is.null(ctx$tg_sens)) qsp2_cal_I(st, ctx$tg_sens) else NULL
      if (!is.null(Is)) flags <- cbind(flags, with_C5_grade = r13 & passes(cbind(I, Is)))
      Irec <- if (is.null(Is)) I else cbind(I, Is)  # recorded with the I_ prefix; not in any rule below
      c2 <- c("C2_p_combo", "C2_p_mono", "C2_ratio")
      I3 <- I
      st3 <- st[, c("C2T3_p_combo", "C2T3_p_mono", "C2T3_ratio"), drop = FALSE]
      colnames(st3) <- c2
      I3[, c2] <- qsp2_cal_I(st3, ctx$tg[match(c2, ctx$tg$id), ])
      flags <- cbind(flags, C2_as_T3 = r13 & passes(I3))
      c1 <- ids[srcs == "C1"]
      I1 <- I
      st1 <- st[, sub("^C1_", "C1mix_", c1), drop = FALSE]
      colnames(st1) <- c1
      I1[, c1] <- qsp2_cal_I(st1, ctx$tg[match(c1, ctx$tg$id), ])
      flags <- cbind(flags, C1_mixture = r13 & passes(I1))
      colnames(flags) <- paste0("pass_", colnames(flags))
      keep <- which(rowSums(flags) > 0)
      b <- which.min(ifelse(lost_ok, t2a[, 1], Inf))
      out$best <- c(tri_row = b, unlist(ctx$tri[b, ]), Imax = unname(t2a[b, 1]), I2nd = unname(t2a[b, 2]),
                    R1 = qsp2_req_R1(st)[b], R3 = qsp2_req_R3(st)[b], st[b, ],
                    setNames(Irec[b, ], paste0("I_", colnames(Irec))))
      out$keep <- if (length(keep)) cbind(tri_row = keep, as.matrix(ctx$tri[keep, ]),
                                          Imax = unname(t2a[keep, 1]), I2nd = unname(t2a[keep, 2]),
                                          flags[keep, , drop = FALSE] * 1,
                                          st[keep, , drop = FALSE],
                                          setNames(as.data.frame(Irec[keep, , drop = FALSE]),
                                                   paste0("I_", colnames(Irec))) |> as.matrix()) else NULL
    }
    out
  }, error = function(e) list(id = job$id, status = paste("error:", conditionMessage(e))))
  res
}

#' C6's troponin T statistic for one member at several regimen mixes, because Pereyra
#' Pietri does not report its mix and 10% is assumed.  job$draw holds the parameters and
#' job$tri the member's population.  The statistic is each patient's highest troponin, the
#' same quantity the calibration fits.
qsp2_cal_tropmix <- function(job, ctx) {
  tryCatch({
    pars <- qsp2_cal_pars(job$draw)
    run <- qsp2_cal_run(pars, qsp2_cal_cohorts()["S0T"], ctx$sgrid)
    WM <- qsp2_pop_W(run$tr$mono, job$tri); WC <- qsp2_pop_W(run$tr$combo, job$tri)
    hm <- WM %*% run$obs$S0T$mono$troppk; hc <- WC %*% run$obs$S0T$combo$troppk
    vapply(ctx$mix, function(w) qsp2_pop_tropmed_hist((1 - w) * hm + w * hc)[1], 0)
  }, error = function(e) rep(NA_real_, length(ctx$mix)))
}
