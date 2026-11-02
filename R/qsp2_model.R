# ---------------------------------------------------------------------------
# qsp2_model.R  -  the model: parameters, dosing calendars, and the simulator.
#
# Implements docs/model_specification.md.  The
# right-hand side lives in src/qsp2_rhs.c; this file builds its parameter vector,
# sets up the states, and integrates piecewise between dose boundaries so that
# each infusion is delivered in full.
#
# Definitions only: sourcing this file runs nothing.
# ---------------------------------------------------------------------------

# --- fixed parameters, every one from data/extraction.csv --------------------------
# Pharmacokinetics.  The sources report nivolumab and ipilimumab clearances in mL/h
# and L/h; they are converted to L/day here and block B4 checks each conversion.
QSP2_PK <- list(
  nivolumab     = list(CL = 9.4 * 24 / 1000, Vc = 3.63, Q = 32.1 * 24 / 1000, Vp = 2.78,
                       aCL = 0.566, aV = 0.597, BWref = 80,   src = "bajaj2017"),
  ipilimumab    = list(CL = 0.0150 * 24,     Vc = 4.15, Q = 0.0411 * 24,      Vp = 3.11,
                       aCL = 0.642, aV = 0.708, BWref = 80,   src = "feng2014"),
  pembrolizumab = list(CL = 0.22,            Vc = 3.48, Q = 0.795,            Vp = 4.06,
                       aCL = 0.595, aV = 0.489, BWref = 76.8, src = "ahamadi2017")
)

# Troponin, from Kristensen 2024: half-lives of 134.1 min (T) and 239.7 min (I).
QSP2_TROPONIN <- list(
  hs_cTnT = list(k_el = log(2) / (134.1 / 1440), limit = 14.0,   Vd = 5.5, cv = 0.09),
  hs_cTnI = list(k_el = log(2) / (239.7 / 1440), limit = 26.2,   Vd = 5.5, cv = 0.09)
)

QSP2_MYO <- list(Mlv = 124, Vlv_mL = 118, EF0 = 60, EFmin = 10)

# Tier T2, from the specification Section 9: the biopsy threshold of 7 CD3+ cells per
# mm2, converted with the Abercrombie correction (section 4 um, nuclear diameter 7 um).
# This gives the whole-ventricle equivalent of the density Caforio 2013 proposes for an
# endomyocardial biopsy specimen.
QSP2_T2 <- list(cells_per_mm2 = 7, t_sec_um = 4, d_nuc_um = 7)
QSP2_T2$cells <- QSP2_T2$cells_per_mm2 /
  ((QSP2_T2$t_sec_um + QSP2_T2$d_nuc_um) / 1000) * QSP2_MYO$Vlv_mL * 1000

#' The infiltrate at which the model calls a patient a case of myocarditis.
#'
#' The cohorts this model is compared with diagnose myocarditis by consensus clinical
#' criteria, not by biopsy, and the biopsy density is read in a sampled fragment rather
#' than averaged over the ventricle.  The tier is therefore `f_dx` times the whole-organ
#' biopsy equivalent, with `f_dx` fitted on calibration data, because the conversion is
#' uncertain (Section 9 of the specification).  `f_dx` may not be read as a focality fraction: it absorbs the
#' arithmetic of the conversion, the focality of the infiltrate and the gap between a
#' histological threshold and a clinical diagnosis, and those are not separately
#' identified.
qsp2_t2_cells <- function(pars) {
  f <- if ("f_dx" %in% names(pars)) pars[["f_dx"]] else 1
  QSP2_T2$cells * f
}

#' Ejection fraction from the tissue state and the inflammatory burden (Section 7 of the
#' specification).  `sat` is the saturating burden eff / (Kh + eff), 0 for tissue alone.
qsp2_lvef <- function(v, d, pars, EF0 = QSP2_MYO$EF0, sat = 0) QSP2_MYO$EFmin +
  (EF0 - QSP2_MYO$EFmin) * (v + pars[["eta"]] * d)^pars[["pexp"]] *
  (1 - qsp2_phi_inf(pars) * sat)

#' The greatest fraction of function inflammation can take away.  Absent from a parameter
#' vector means zero.
qsp2_phi_inf <- function(pars) if ("phi_inf" %in% names(pars)) pars[["phi_inf"]] else 0

#' The saturating cardiac effector burden eff / (Kh + eff) of a simulation result.  `eff` is the
#' burden that drives killing in module M4, so the depression of function saturates where killing
#' saturates and no second half-saturation constant is introduced.
qsp2_inflam_sat <- function(r, pars) {
  eff <- r$Eh + pars[["wL"]] * r$L + pars[["wM"]] * r$Mq
  # The depression saturates where killing saturates, at Kh, so that it adds one fitted
  # parameter rather than two.  A vector carrying K_inf uses that instead (an alternative, not
  # used in the calibration).
  k <- if ("K_inf" %in% names(pars) && is.finite(pars[["K_inf"]])) pars[["K_inf"]] else pars[["Kh"]]
  eff / (k + eff)
}

# --- the full parameter vector -------------------------------------------------------
# The order here is the order the C code expects; changing one without the other is
# the classic silent error, so B4 checks the length and a term-by-term agreement
# between the C right-hand side and an independent one written in R.
QSP2_PAR_NAMES <- c(
  "CLn", "Vcn", "Qn", "Vpn", "CLi", "Vci", "Qi", "Vpi", "CLp", "Vcp", "Qp", "Vpp",
  "EC50pd1", "EC50ctla4", "lam_p", "gC", "k_act", "del_P", "rho", "gP", "Emax", "del_E",
  "k_tr", "phi", "pi_h", "gF", "Kg", "del_h", "eps_x", "k_lic", "K_lic", "del_L",
  "k_m", "Mmax", "del_m", "lam_g", "k_g", "del_g", "k_x", "del_x", "k_kill", "Kh", "wL", "NT",
  "wM", "Sus", "r_rep", "r_scar", "theta", "fc", "k_rel", "k_el", "Vd", "Mlv", "Tn0",
  "infn", "infi", "infp", "gPp", "Svul")

#' Default parameters: fixed values from the sources, fitted ones at the middle of
#' their prior. These are a starting point for verification, not a calibration.
qsp2_default_pars <- function(BW = 80, assay = "hs_cTnT") {
  pk <- QSP2_PK; tr <- QSP2_TROPONIN[[assay]]
  sc <- function(d, what) {
    r <- pk[[d]]
    if (what == "CL") r$CL * (BW / r$BWref)^r$aCL
    else if (what == "Q") r$Q * (BW / r$BWref)^r$aCL
    else if (what == "Vc") r$Vc * (BW / r$BWref)^r$aV
    else r$Vp * (BW / r$BWref)^r$aV
  }
  c(CLn = sc("nivolumab", "CL"), Vcn = sc("nivolumab", "Vc"),
    Qn = sc("nivolumab", "Q"),   Vpn = sc("nivolumab", "Vp"),
    CLi = sc("ipilimumab", "CL"), Vci = sc("ipilimumab", "Vc"),
    Qi = sc("ipilimumab", "Q"),   Vpi = sc("ipilimumab", "Vp"),
    CLp = sc("pembrolizumab", "CL"), Vcp = sc("pembrolizumab", "Vc"),
    Qp = sc("pembrolizumab", "Q"),   Vpp = sc("pembrolizumab", "Vp"),
    EC50pd1 = 0.5, EC50ctla4 = 1.5,
    lam_p = 1e4, gC = 3, k_act = 0.2, del_P = 0.05,
    rho = 0.25, gP = 3, Emax = 1e9, del_E = 0.05,
    k_tr = 0.25, phi = 0.3, pi_h = 0.05, gF = 2, Kg = 1e5, del_h = 0.1, eps_x = 0.05,
    k_lic = 0.3, K_lic = 1e5, del_L = 0.05,
    k_m = 0.02, Mmax = 1e8, del_m = 0.1,
    lam_g = 1e3, k_g = 0.01, del_g = 0.05,
    k_x = 20, del_x = 0.05,
    k_kill = 1e-4, Kh = 7.51e7, wL = 1, NT = 4,
    wM = 0.2, Sus = 1, r_rep = 0.05, r_scar = 0.02,
    theta = 1.07e5, fc = 0.2, k_rel = 0.1,
    k_el = tr$k_el, Vd = tr$Vd, Mlv = QSP2_MYO$Mlv, Tn0 = 5,
    infn = 0, infi = 0, infp = 0, gPp = 1, Svul = 1,
    # the function map (Section 7) and the symptom rule (Section 9);
    # these are used by the R layer and are not passed to the solver
    eta = 0.5, pexp = 1.5, k_sym = 5, p_bg = 0.005,
    # where the myocarditis tier sits, as a fraction of the whole-ventricle biopsy equivalent
    f_dx = 1,
    # patient variation and the observation model.  None of these reaches the solver.
    #   sd_tau  log-scale spread of the patient's tempo
    #   L_tx    days from diagnosis to the halt of myocyte killing
    #   s_ksym  log-scale spread of the symptom threshold between patients
    #   s_b     log-scale spread of the baseline troponin between patients
    #   lam_bg  rate of background troponin episodes per day
    #   m_bg    median peak of a background episode, ng/L on the hs-cTnT scale
    #   s_bg    log-scale spread of that peak, held at 1.0
    sd_tau = 0.5, L_tx = 3, s_ksym = 0.5, s_b = 1, lam_bg = 0.01, m_bg = 40, s_bg = 1,
    # the greatest fraction of function that inflammation can take away, at the
    # middle of its uniform prior from 0 to 0.8
    phi_inf = 0.4,
    # how strongly the rate of background rises follows the patient's own baseline
    g_bg = 0)
}

# --- tempo ------------------------------------------------------------------------------
# A patient's tempo multiplies the six rates of module M3 together, so the systemic response
# runs at that pace while its size is still set by the precursor supply.  The constraint
# rho < del_E + k_tr is unaffected because both sides scale together.
# The rates at which the heart handles the same response scale with tempo too, because they
# describe the clone inside the heart rather than the tissue.  The myocardium's own r_rep and r_scar, the killing constant
# k_kill and every capacity (K_lic, Kg, Mmax, Kh) keep their own values.
QSP2_TAU_RATES <- c("lam_p", "k_act", "del_P", "rho", "del_E", "k_tr",
                    "pi_h", "del_h", "eps_x", "k_lic", "del_L", "k_m", "del_m",
                    "lam_g", "k_g", "del_g", "k_x", "del_x")
# The grid must be able to hold the sd_tau prior it is given.  The end cells of
# qsp2_pop_cellmass run to minus and plus infinity, so a grid narrower than the distribution
# piles both tails onto its two edge nodes (a sixteen-fold grid would put 41.9% of the
# population on two tempo values at sd_tau = 1.5).  This span keeps the end cells under 1% at
# every sd_tau the prior allows, at a spacing of a factor of sqrt(2).
QSP2_TAU_GRID  <- 2^seq(-6, 6, by = 0.5)          # 25 values, a 4096-fold range

#' A patient's parameter vector at tempo `tau`.
qsp2_apply_tempo <- function(pars, tau) {
  for (nm in QSP2_TAU_RATES) pars[[nm]] <- pars[[nm]] * tau
  pars
}

# --- cohorts -----------------------------------------------------------------------------
# A cohort is its testing calendar together with the assay and threshold that cohort used.
# Predicting every cohort at one threshold would change the rate of raised results by up to
# a factor of two, in both directions.  Each threshold below is a verified row of
# data/extraction.csv.
#
#' Build a cohort description.
#' @param days  the testing calendar, in days from the first dose
#' @param assay "hs_cTnI" or "hs_cTnT"; sets the clearance rate as well as the limit
#' @param limit the cohort's own threshold in ng/L; the assay's own limit when NULL
#' @param note  free text, used where a threshold is not reported by the source
qsp2_cohort <- function(days, assay = "hs_cTnI", limit = NULL, note = "") {
  stopifnot(assay %in% names(QSP2_TROPONIN))
  list(days = sort(unique(days[is.finite(days) & days >= 0])), assay = assay,
       limit = if (is.null(limit)) QSP2_TROPONIN[[assay]]$limit else limit, note = note)
}

#' Initial state. Rm0 is the resident memory pool; everything downstream starts empty.
qsp2_init_state <- function(pars, Rm0 = 1e5, Gr0 = NULL) {
  NT <- as.integer(round(pars[["NT"]]))
  if (is.null(Gr0)) Gr0 <- pars[["lam_g"]] / pars[["del_g"]]   # steady state with no effectors
  # With no drug there is no priming, so the whole chain starts empty and
  # the untreated model sits exactly at rest.
  y <- c(A1n = 0, A2n = 0, A1i = 0, A2i = 0, A1p = 0, A2p = 0,
         P = 0,
         E = 0,
         setNames(rep(0, NT), paste0("T", seq_len(NT))),
         Eh = 0, Rm = Rm0, L = 0, Mq = 0, Gr = Gr0, X = 0,
         v = 1, d = 0, s = 0, Qs = 0, Tn = pars[["Tn0"]])
  y
}

QSP2_NOUT <- 8L
QSP2_OUTNAMES <- c("Cn", "Ci", "Cp", "Opd1", "Octla4", "kappa", "Tcells_heart", "tissue_total")

# --- dosing ---------------------------------------------------------------------------
#' Build a dosing calendar.
#' Each row is one infusion: day, drug, dose in mg, and infusion duration in days.
qsp2_regimen <- function(kind = c("nivo_q2w", "nivo_q4w", "pembro_q3w", "pembro_q6w",
                                  "combo_nivo1_ipi3", "combo_nivo3_ipi1", "ipi_alone", "none"),
                         BW = 80, horizon = 365) {
  kind <- match.arg(kind)
  inf_short <- 0.5 / 24; inf_long <- 1.5 / 24
  mk <- function(days, drug, dose, dur) if (!length(days)) NULL else
    data.table(day = days, drug = drug, dose = dose, dur = dur)
  # a horizon shorter than the start of a phase gives that phase no doses at all; seq() would
  # otherwise count backwards and fail (found by a block B4 test of a finite induction course)
  sq <- function(from, to, by) if (to < from) numeric(0) else seq(from, to, by = by)
  out <- switch(kind,
    none = NULL,
    nivo_q2w   = mk(sq(0, horizon, 14), "nivolumab", 240, inf_short),
    nivo_q4w   = mk(sq(0, horizon, 28), "nivolumab", 480, inf_short),
    pembro_q3w = mk(sq(0, horizon, 21), "pembrolizumab", 200, inf_short),
    pembro_q6w = mk(sq(0, horizon, 42), "pembrolizumab", 400, inf_short),
    ipi_alone  = mk(sq(0, min(63, horizon), 21), "ipilimumab", 3 * BW, inf_long),
    combo_nivo1_ipi3 = rbind(
      mk(sq(0, min(63, horizon), 21), "nivolumab", 1 * BW, inf_short),
      mk(sq(0, min(63, horizon), 21), "ipilimumab", 3 * BW, inf_long),
      mk(sq(84, horizon, 14), "nivolumab", 240, inf_short)),
    combo_nivo3_ipi1 = rbind(
      mk(sq(0, min(63, horizon), 21), "nivolumab", 3 * BW, inf_short),
      mk(sq(0, min(63, horizon), 21), "ipilimumab", 1 * BW, inf_long),
      mk(sq(84, horizon, 14), "nivolumab", 240, inf_short))
  )
  if (is.null(out)) return(data.table(day = numeric(0), drug = character(0),
                                      dose = numeric(0), dur = numeric(0)))
  out[day <= horizon][order(day)]
}

# --- the simulator ----------------------------------------------------------------------
#' Integrate the model, piecewise between dose boundaries.
#'
#' Infusions are zero-order inputs, so within any interval where the set of running
#' infusions does not change the right-hand side is smooth.  The interval boundaries
#' are every infusion start and end; the solver is restarted at each, which is what
#' delivers every dose in full.
qsp2_simulate <- function(pars, regimen, times, y0 = NULL, assay = NULL,
                          rtol = 1e-10, atol = 1e-12) {
  stopifnot(all(QSP2_PAR_NAMES %in% names(pars)))
  # `assay` sets the troponin clearance rate and volume.  Troponin T clears about twice as
  # fast as troponin I, so the two assays do not see the same transient.
  #
  # The default is NULL, meaning "use the rates in `pars`".  It has to be NULL rather
  # than an assay name, because the B4 dimensional check runs the whole model in hours by
  # dividing every rate by 24, and an assay named here would silently put `k_el` back
  # into days and make that check fail.  Naming an assay is therefore a deliberate act.
  if (!is.null(assay)) {
    stopifnot(assay %in% names(QSP2_TROPONIN))
    pars[["k_el"]] <- QSP2_TROPONIN[[assay]]$k_el
    pars[["Vd"]]   <- QSP2_TROPONIN[[assay]]$Vd
  }
  if (is.null(y0)) y0 <- qsp2_init_state(pars)
  load_qsp2_dll()

  ev <- if (nrow(regimen)) rbind(
      data.table(t = regimen$day, drug = regimen$drug, rate = regimen$dose / regimen$dur),
      data.table(t = regimen$day + regimen$dur, drug = regimen$drug,
                 rate = -regimen$dose / regimen$dur)) else
      data.table(t = numeric(0), drug = character(0), rate = numeric(0))

  bounds <- sort(unique(c(min(times), max(times), ev$t, times)))
  bounds <- bounds[bounds >= min(times) & bounds <= max(times)]
  segs <- sort(unique(c(min(times), ev$t[ev$t > min(times) & ev$t < max(times)], max(times))))

  rate <- c(nivolumab = 0, ipilimumab = 0, pembrolizumab = 0)
  res <- NULL
  y <- y0
  for (k in seq_len(length(segs) - 1L)) {
    t0 <- segs[k]; t1 <- segs[k + 1L]
    now <- ev[abs(t - t0) < 1e-12]
    if (nrow(now)) for (i in seq_len(nrow(now)))
      rate[[now$drug[i]]] <- rate[[now$drug[i]]] + now$rate[i]
    rate[abs(rate) < 1e-9] <- 0                      # kill floating-point residue

    pv <- pars
    pv[["infn"]] <- rate[["nivolumab"]]
    pv[["infi"]] <- rate[["ipilimumab"]]
    pv[["infp"]] <- rate[["pembrolizumab"]]

    tt <- sort(unique(c(t0, times[times > t0 & times < t1], t1)))
    out <- deSolve::lsoda(y = y, times = tt, func = "qsp2_derivs", parms = pv[QSP2_PAR_NAMES],
                          dllname = "qsp2_rhs", initfunc = "qsp2_init",
                          nout = QSP2_NOUT, outnames = QSP2_OUTNAMES,
                          rtol = rtol, atol = atol)
    y <- out[nrow(out), -1][seq_along(y0)]
    names(y) <- names(y0)
    res <- rbind(res, if (is.null(res)) out else out[-1, , drop = FALSE])
  }
  res <- as.data.table(res)
  setnames(res, "time", "t")
  res[t %in% times]
}


# --- surveillance schedules --------------------------------------------------------------
#' The days on which troponin is measured, for a schedule and a regimen.
#'
#' Schedules tied to doses are converted to calendar days using that regimen's own
#' dosing calendar, so "before dose 2" is day 14 on a two-weekly regimen and day 21 on
#' a three-weekly one.
qsp2_schedule <- function(which = c("S0", "S1", "S2", "S3", "S4"), regimen, horizon = 365) {
  which <- match.arg(which)
  dose_days <- sort(unique(regimen$day))
  if (!length(dose_days)) dose_days <- 0
  days <- switch(which,
    # baseline only; anything later is found by symptoms, not by testing
    S0 = 0,
    # 2022 ESC guideline: baseline, before doses 2, 3 and 4, then every third dose
    S1 = {
      d <- c(0, dose_days[2:min(4, length(dose_days))])
      # every third dose from the fifth; a calendar with fewer than five doses has none
      rest <- if (length(dose_days) >= 5) dose_days[seq(5, length(dose_days), by = 3)] else numeric(0)
      sort(unique(c(d, rest)))
    },
    # before every infusion
    S2 = sort(unique(c(0, dose_days))),
    # weekly for the first eight weeks, then as the guideline
    S3 = {
      g <- qsp2_schedule("S1", regimen, horizon)
      sort(unique(c(seq(0, 56, by = 7), g[g > 56])))
    },
    # every four weeks for six months, then every three months
    # (a horizon shorter than 252 days has no quarterly tests)
    S4 = sort(unique(c(seq(0, 168, by = 28), if (horizon >= 252) seq(252, horizon, by = 91)))))
  days <- days[is.finite(days) & days >= 0 & days <= horizon]
  sort(unique(days))
}
