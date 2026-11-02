# ---------------------------------------------------------------------------
# lib_pop.R  -  the virtual population engine.
#
# The model's equations are solved on a grid of the precursor supply S crossed with the
# patient's tempo tau.  Everything else that differs between patients, and everything about
# how a cohort measures them, acts on a stored trajectory without solving it again:
#   * the symptom threshold k_sym_i.  It matters only through the day a patient presents with
#     troponin, whose distribution at each grid point is a normal probability of the running
#     maximum of troponin;
#   * the baseline LVEF EF0_i.  It matters only through the value of G at which the patient
#     meets the LVEF criterion (below 50% and at least 10 points below baseline), so the day
#     of presentation by LVEF, T3 given a nadir, the nadir-LVEF bands and the LVEF arm of the
#     C3 definition are all probabilities of intervals of a normal EF0, in closed form;
#   * the two presentation routes compete.  Their days are cut into equal-probability slices
#     placed at their median days (eight for troponin, two for LVEF), plus the exact
#     probability of not presenting;
#   * the baseline troponin b_i (16 points, thin ones on either side of the cohort's limit
#     where measurement error decides; the share with a raised baseline is exact);
#   * the measurement error of the assay and background troponin episodes (closed form);
#   * diagnosis, treatment after diagnosis, and severity over the treated course.
# Population shares are weighted sums over the grid.  No outcome is assumed to be monotone.
#
# Why not fixed nodes.  Most outcomes here are thresholds.  Gauss-Hermite nodes
# are exact for smooth integrands and biased for thresholds (the lowest 3-point node for LVEF,
# 47.9%, carries a sixth of the weight and sits below the 50% cut-off), and fixed nodes for the
# symptom threshold miss the late presenters in its tail, who are the ones that become severe.
# A Monte Carlo of individual patients checks the engine (block B4).
#
# LVEF criteria.  With baselines spread around 60%, "LVEF below 50%" alone would
# make a patient starting at 51% severe after a one-point fall, and a patient starting below 50%
# severe with no fall at all.  T3, and presentation by LVEF, therefore need a fall of at least
# 10 points to below 50%, the LVEF-only criterion of moderate cardiac dysfunction in the 2022
# ESC guideline (Table 3).  For a patient starting at 60% nothing changes.
#
# Definitions only: sourcing this file runs nothing.
# ---------------------------------------------------------------------------

QSP2_REF_URL   <- 26.2                              # hs-cTnI (Abbott): the reference assay
QSP2_TNT_URL   <- 14.0                              # hs-cTnT (Roche)
# Clearance takes hours and release changes over days, so troponin T follows troponin I
# in proportion to the ratio of clearance rates; block B4 measures the error.
QSP2_TNT_RATIO <- QSP2_TROPONIN$hs_cTnI$k_el / QSP2_TROPONIN$hs_cTnT$k_el
QSP2_CV_ANA    <- 0.09                              # assay coefficient of variation (C6 source)
QSP2_SDLOG     <- sqrt(log(1 + QSP2_CV_ANA^2))
QSP2_WORKUP    <- 7                                 # days allowed for the diagnostic work-up
QSP2_EF_CUT    <- 50                                # T3 and the ejection-fraction symptom rule
QSP2_EF_FALL   <- 10                                # points of fall that T3 and presentation need
QSP2_EF_SD     <- 7                                 # EF0 normal, mean 60%, SD 7%
QSP2_C3_EF     <- 53                                # C3 (Nishikawa 2022): a fall of >10 to <53%
QSP2_TROP_BINS <- seq(0, 5, by = 0.05)              # log10 ng/L, for troponin at diagnosis
QSP2_TP_SLICES <- 8L                                # slices of the troponin-presentation days
QSP2_TE_SLICES <- 2L                                # slices of the LVEF-presentation days

# The LVEF criterion at G = g holds for baselines e in the interval I(g) = (L(g), U(g)):
# a fall of at least 10 points needs e > EFmin + 10 / (1 - g), and a value below 50% needs
# e < EFmin + 40 / g.  I(g) shrinks as g rises and is empty from G_MAX = 0.8, reached at e = 60.
QSP2_G_MAX <- (QSP2_EF_CUT - QSP2_MYO$EFmin) / (QSP2_EF_CUT - QSP2_MYO$EFmin + QSP2_EF_FALL)
QSP2_E_MID <- QSP2_EF_CUT + QSP2_EF_FALL

#' The value of G below which a patient with baseline LVEF e meets the LVEF criterion.
qsp2_g_crit <- function(e) {
  em <- QSP2_MYO$EFmin
  pmax(pmin((QSP2_EF_CUT - em) / (e - em), 1 - QSP2_EF_FALL / (e - em)), 0)
}
.efL <- function(g) ifelse(g < QSP2_G_MAX,
                           QSP2_MYO$EFmin + QSP2_EF_FALL / (1 - pmin(g, QSP2_G_MAX)), QSP2_E_MID)
.efU <- function(g) ifelse(g < QSP2_G_MAX,
                           QSP2_MYO$EFmin + (QSP2_EF_CUT - QSP2_MYO$EFmin) / pmax(g, 1e-300),
                           QSP2_E_MID)
#' P(a < EF0 < b) for EF0 normal(60, 7), elementwise.
.pe <- function(a, b) pmax(pnorm((b - QSP2_MYO$EF0) / QSP2_EF_SD) -
                             pnorm((a - QSP2_MYO$EF0) / QSP2_EF_SD), 0)
#' P(the patient's critical value exceeds g) = P(EF0 in I(g)): the share that meets the LVEF
#' criterion once G has fallen to g.
qsp2_S_c <- function(g) .pe(.efL(g), .efU(g))
#' P(EF0 in [p, q) and critical value in (lo, hi]), i.e. EF0 in [p, q) and in I(lo) \ I(hi).
.p_in_set <- function(p, q, lo, hi)
  .pe(pmax(p, .efL(lo)), pmin(q, .efL(hi))) + .pe(pmax(p, .efU(hi)), pmin(q, .efU(lo)))

#' Severity given the nadir of G, for patients whose critical value lies in (lo, hi]: T3 and
#' the three nadir-LVEF bands (40 to 49%, 30 to 39%, below 30%), as shares of that set.
#' `forced` marks patients already known to have met the criterion (they presented by LVEF),
#' who are T3 whatever the nadir; their band is read from the nadir, and a nadir that would
#' leave them at or above 50% puts them in the mildest band.  Returns a matrix with 4 columns.
.sev_c <- function(gs, lo, hi, forced = FALSE) {
  n <- max(length(gs), length(lo), length(hi))
  gs <- rep_len(pmin(pmax(gs, 1e-12), 1), n); lo <- rep_len(lo, n); hi <- rep_len(hi, n)
  forced <- rep_len(forced, n)
  em <- QSP2_MYO$EFmin
  pset <- qsp2_S_c(lo) - qsp2_S_c(hi)
  m <- pmax(gs, lo)
  in3 <- m < hi
  b40c <- em + (40 - em) / gs; b30c <- em + (30 - em) / gs
  inT3 <- function(p, q) .pe(pmax(p, .efL(m)), pmin(q, .efL(hi))) +
    .pe(pmax(p, .efU(hi)), pmin(q, .efU(m)))
  t3  <- ifelse(forced, pset, ifelse(in3, qsp2_S_c(m) - qsp2_S_c(hi), 0))
  b30 <- ifelse(forced, .p_in_set(b30c, b40c, lo, hi), ifelse(in3, inT3(b30c, b40c), 0))
  b0  <- ifelse(forced, .p_in_set(-Inf, b30c, lo, hi), ifelse(in3, inT3(-Inf, b30c), 0))
  out <- cbind(t3 = t3, b40 = t3 - b30 - b0, b30 = b30, b0 = b0) / ifelse(pset > 0, pset, 1)
  out[!(pset > 0), ] <- 0
  pmin(pmax(out, 0), 1)
}

.first_day <- function(z) { w <- which(z); if (length(w)) w[1] - 1L else NA_integer_ }

# --- trajectories on the (S, tau) grid ------------------------------------------------------

#' Solve the model at every grid point and store what the observation layer needs.
#'
#' @return a list: `nodes` (S, tau and their grid indices), daily matrices `X` (troponin I
#'   released by drug-induced injury, ng/L, i.e. Tn - Tn0), `V`, `D` and `G`
#'   (= (v + eta d)^pexp, so LVEF = EFmin + (EF0 - EFmin) G), `cumG` (running minimum of G),
#'   and for each node `t2c` (first day of T2) and `t1ref` (first day of T1 on the reference
#'   assay with the reference baseline, which is what requirement R1 reads).
qsp2_pop_traj <- function(pars, regimen, sgrid, taugrid = QSP2_TAU_GRID, H = 365) {
  days <- 0:H
  nodes <- expand.grid(i_s = seq_along(sgrid), i_t = seq_along(taugrid))
  nodes$S <- sgrid[nodes$i_s]; nodes$tau <- taugrid[nodes$i_t]
  N <- nrow(nodes)
  X <- V <- D <- CELLS <- SAT <- matrix(NA_real_, N, H + 1L)
  ok <- logical(N)
  for (k in seq_len(N)) {
    p <- qsp2_apply_tempo(pars, nodes$tau[k]); p[["Sus"]] <- nodes$S[k]
    r <- tryCatch(qsp2_simulate(p, regimen, days, assay = "hs_cTnI"), error = function(e) NULL)
    # a truncated solve would drop rows and silently shift every test day, so it is refused
    if (is.null(r) || nrow(r) != H + 1L || anyNA(r$Tn) || anyNA(r$v)) next
    X[k, ] <- pmax(r$Tn - p[["Tn0"]], 0)
    V[k, ] <- r$v; D[k, ] <- r$d; CELLS[k, ] <- r$Tcells_heart
    SAT[k, ] <- qsp2_inflam_sat(r, p)
    ok[k] <- TRUE
  }
  # function falls with lost tissue and, reversibly, with the inflammatory burden.
  # Gtis is the tissue-only value the course returns to once the infiltrate clears.
  Gtis <- (V + pars[["eta"]] * D)^pars[["pexp"]]
  G <- Gtis * (1 - qsp2_phi_inf(pars) * SAT)
  cumG <- G
  if (any(ok)) cumG[ok, ] <- t(apply(G[ok, , drop = FALSE], 1, cummin))
  tier <- qsp2_t2_cells(pars)
  # T2 is T1 together with the infiltrate above the tier; both are read on the reference assay
  # with the reference baseline, on the same day
  t1on  <- X + pars[["Tn0"]] > QSP2_REF_URL & !is.na(X)
  t2c   <- apply(CELLS > tier & !is.na(CELLS) & t1on, 1, .first_day)
  t1ref <- apply(t1on, 1, .first_day)
  list(pars = pars, H = H, sgrid = sgrid, taugrid = taugrid, nodes = nodes,
       X = X, V = V, D = D, SAT = SAT, G = G, Gtis = Gtis, cumG = cumG, t2c = t2c,
       t1ref = t1ref, ok = ok, peakX = apply(X, 1, max), regimen = regimen)
}

#' n equal-probability slices of a mass over days 0..H: each slice's median day, first and
#' last day, and the total.
.slices <- function(mass, n) {
  P <- sum(mass)
  if (!(P > 0)) return(list(rep = rep(NA_integer_, n), a = rep(NA_integer_, n),
                            b = rep(NA_integer_, n), P = 0))
  cm <- cumsum(mass)
  at <- function(q) findInterval(q, cm, left.open = TRUE)   # first day the cumulative reaches q
  i <- seq_len(n)
  list(rep = at(P * (i - 0.5) / n), a = at(P * (i - 1) / n + P * 1e-9),
       b = at(pmin(P * i / n, P * (1 - 1e-9))), P = P)
}

#' What does not depend on the cohort, for a follow-up of H days.
#'
#' Patient classes at each node: presentation by troponin (n_t slices), by LVEF (n_e slices),
#' and none within H.  For each class: its median presentation day (TP), its probability (WM),
#' and the set (LO, HI] its patients' critical values lie in.
#'  * Troponin: presented by day t with probability Phi((log(M(t) / URL) - log k_sym) / s_ksym),
#'    where M(t) is the running maximum of the drug-released troponin on top of the reference
#'    baseline Tn0; 0 while M(t) is at or below the limit.  A chronically raised baseline is not
#'    a symptom, so the patient's own baseline does not enter.  Their LVEF has not met the
#'    criterion when they present, so their critical value is at most cumG at that day.
#'  * LVEF: the criterion is met on day t by patients whose critical value lies in
#'    (cumG(t), cumG(t - 1)], and they have not presented by troponin before day t.
#'  * None: neither by day H.
#' Also stored: the nadir of G within the follow-up if killing halts on each day (gpost).
#' @param persist  an alternative assumption, not used in the calibration: after the halt of
#'   killing the inflammatory depression is not taken to resolve with the infiltrate but follows the course
#'   the infiltrate would have taken anyway, so the nadir is searched over the whole follow-up.
qsp2_pop_prep <- function(traj, pp, H = traj$H, n_t = QSP2_TP_SLICES, n_e = QSP2_TE_SLICES,
                          persist = FALSE) {
  pp <- as.list(pp)
  stopifnot(H <= traj$H)
  N <- nrow(traj$nodes)
  cls <- c(rep("T", n_t), rep("E", n_e), "N"); nM <- length(cls)
  # After a halt on day t, v + eta d
  # moves monotonically from its value on day t towards v + d r_rep / (r_rep + r_scar), and on
  # day H it is v_t + d_t (rho + (eta - rho) exp(-(r_rep + r_scar) (H - t))).  The nadir within
  # the follow-up is the smaller of the running minimum up to t and that value; a halt after
  # day H leaves the untreated nadir up to H.  Measured this way, a treated course is never
  # worse than an untreated one and an earlier halt never worse than a later one, at every day.
  rr <- pp$r_rep + pp$r_scar; rho_rep <- pp$r_rep / rr
  eta <- traj$pars[["eta"]]; pexp <- traj$pars[["pexp"]]
  tt <- 0:traj$H
  frac <- rho_rep + (eta - rho_rep) * exp(-rr * pmax(H - tt, 0))
  # After the halt the immunosuppression that stops killing also resolves the infiltrate, so
  # the value the course returns to is computed from tissue alone.
  gpost <- pmin(traj$cumG, (traj$V + sweep(traj$D, 2, frac, "*"))^pexp)
  if (persist) {
    # the depression follows the infiltrate, which killing does not drive, so the stored burden
    # is still the right one after the halt; the nadir is the smallest value the course reaches
    phi <- qsp2_phi_inf(traj$pars)
    for (k in which(traj$ok)) {
      gp <- numeric(traj$H + 1L)
      for (th in 0:H) {
        s <- th:H
        vv <- traj$V[k, th + 1L] + traj$D[k, th + 1L] *
          (rho_rep + (eta - rho_rep) * exp(-rr * (s - th)))
        gp[th + 1L] <- min(traj$cumG[k, th + 1L], min(vv^pexp * (1 - phi * traj$SAT[k, s + 1L])))
      }
      if (H < traj$H) gp[(H + 2L):(traj$H + 1L)] <- traj$cumG[k, H + 1L]
      gpost[k, ] <- gp
    }
  }
  if (any(tt > H)) gpost[, tt > H] <- traj$cumG[, H + 1L]
  TP <- matrix(NA_integer_, N, nM); WM <- LO <- HI <- matrix(0, N, nM)
  sk <- max(pp$s_ksym, 1e-4)
  for (k in which(traj$ok)) {
    cg <- traj$cumG[k, seq_len(H + 1L)]
    mx <- cummax(traj$X[k, seq_len(H + 1L)] + pp$Tn0)
    FT <- ifelse(mx > QSP2_REF_URL, pnorm((log(mx / QSP2_REF_URL) - log(pp$k_sym)) / sk), 0)
    SE <- qsp2_S_c(cg)                                   # P(LVEF criterion met by day t)
    mT <- diff(c(0, FT)) * (1 - SE)
    mE <- diff(c(0, SE)) * (1 - c(0, FT[-(H + 1L)]))
    sT <- .slices(mT, n_t); sE <- .slices(mE, n_e)
    cgm <- function(t) ifelse(is.na(t), 1, ifelse(t < 0, 1, cg[pmax(t, 0L) + 1L]))
    TP[k, ] <- c(sT$rep, sE$rep, NA_integer_)
    WM[k, ] <- c(rep(sT$P / n_t, n_t), rep(sE$P / n_e, n_e), (1 - FT[H + 1L]) * (1 - SE[H + 1L]))
    LO[k, ] <- c(rep(0, n_t), ifelse(is.na(sE$b), 0, cgm(sE$b)), 0)
    HI[k, ] <- c(cgm(sT$rep), cgm(sE$a - 1L), cg[H + 1L])
  }
  list(TP = TP, WM = WM, LO = LO, HI = HI, cls = cls, forced = cls == "E", gpost = gpost,
       pp = pp, H = H)
}

# --- cohorts ---------------------------------------------------------------------------------

#' Describe how a cohort measures its patients.
#' @param days     scheduled test days after the baseline, in days from the first dose
#' @param kin      "TnI" or "TnT": the clearance of the assay
#' @param url      the assay's upper reference limit in its own ng/L
#' @param limit    the cohort's threshold for a positive result, in the assay's ng/L
#' @param rule     "default" (a raised baseline needs r_rise times baseline), "double"
#'                 (at least abs_min and at least twice baseline, for everyone),
#'                 "conversion" (a raised baseline cannot convert), "normal_only" (only
#'                 patients with a normal baseline are counted), or "threshold" (every value
#'                 above the limit is positive, whatever the baseline)
#' @param count_baseline  TRUE when the baseline sample itself is counted as a result
#' @param base     NULL, or the cohort's own baseline stratum: list(med, sdlog)
#'                 on the assay's scale, or list(praised) for a reported share above `limit`
#' @param H        follow-up in days, at most the simulated horizon
#' @param w_combo  share of the cohort on combination therapy
#' @param conv     how a troponin I assay other than the reference is related to it
#'                 "url_ratio", the main analysis, puts every
#'                 troponin I concentration (baseline, drug-released, background) on the cohort's
#'                 assay by the ratio of its 99th-percentile upper reference limit to the
#'                 reference assay's, so a result at k times one assay's limit is at k times the
#'                 other's; "none" (an alternative, not used in the calibration) reads ng/L
#'                 as the same on every troponin I assay.  Troponin T cohorts are unaffected.
#' @param expo     the cohort's testing exposure: NULL, every scheduled test is
#'                 done unless a diagnosis comes first; or list(mean, max), the number of tests
#'                 after the baseline is geometric on 1..max with that mean (the patient stops
#'                 treatment, and so testing, after that many doses)
qsp2_pop_cohort <- function(name, days = numeric(0), kin = "TnI", url = QSP2_REF_URL,
                            limit = url, rule = "default", r_rise = 1.5, abs_min = NA_real_,
                            count_baseline = FALSE, base = NULL, H = 365, w_combo = 0,
                            note = "", conv = "url_ratio", expo = NULL) {
  stopifnot(kin %in% c("TnI", "TnT"),
            rule %in% c("default", "double", "conversion", "normal_only", "threshold"),
            rule != "double" || is.finite(abs_min), conv %in% c("url_ratio", "none"),
            is.null(expo) || (expo$mean >= 1 && expo$max >= expo$mean))
  days <- sort(unique(days[is.finite(days) & days > 0 & days <= H]))
  list(name = name, days = days, kin = kin, url = url, limit = limit, rule = rule,
       r_rise = r_rise, abs_min = abs_min, count_baseline = count_baseline, base = base,
       H = H, w_combo = w_combo, note = note, conv = conv, expo = expo)
}

#' A cohort whose patients follow different calendars, e.g. the drug mix of C4:
#' each part is a cohort from qsp2_pop_cohort with the same follow-up and regimen mix, and
#' `w` the share of patients on it.  Shares and day masses are linear in the patients, so the
#' mixture's observation is the weighted sum of the parts' (qsp2_pop_obs_any).
qsp2_pop_cohort_mix <- function(name, parts, w) {
  stopifnot(length(parts) == length(w), abs(sum(w) - 1) < 1e-9, all(w >= 0),
            length(unique(vapply(parts, function(p) p$H, 0))) == 1L,
            length(unique(vapply(parts, function(p) p$w_combo, 0))) == 1L)
  list(name = name, parts = parts, w_parts = w, H = parts[[1]]$H, w_combo = parts[[1]]$w_combo,
       kin = parts[[1]]$kin)
}

#' The factor that puts a reference-assay troponin I concentration on the cohort's troponin I
#' assay.  For troponin T cohorts this is the ratio of the two limits, which is only used for
#' the baseline; drug-released troponin T uses the clearance ratio.
.assay_ratio <- function(co) {
  if (co$kin == "TnT") return(co$url / QSP2_REF_URL)
  if (identical(co$conv, "none")) 1 else co$url / QSP2_REF_URL
}

#' P(K = k), k = 0..nt, for the number of tests after the baseline that a patient's treatment
#' lasts.  NULL exposure puts all the mass on nt.  A geometric on 1..max with the
#' stated mean (the truncation at max is allowed for when solving for its parameter), with the
#' mass from nt onwards placed on nt because only nt tests fall within the follow-up.
qsp2_pop_expo_pk <- function(expo, nt) {
  if (is.null(expo) || nt == 0L) return(c(rep(0, nt), 1))
  mx <- expo$max; k <- seq_len(mx)
  pg <- function(p) { w <- p * (1 - p)^(k - 1); w / sum(w) }
  p <- if (expo$mean >= (mx + 1) / 2 - 1e-9) 1e-9 else
    uniroot(function(p) sum(k * pg(p)) - expo$mean, c(1e-9, 1 - 1e-9), tol = 1e-12)$root
  w <- pg(p)
  out <- numeric(nt + 1L)
  for (i in k) { j <- min(i, nt); out[j + 1L] <- out[j + 1L] + w[i] }
  out
}

#' The baseline points for a cohort, each at the conditional mean of log b within its slice of the
#' baseline's normal z scale; the slices meet at the cohort's limit, so the share with a raised
#' baseline is exact.  The baseline scales between assays by the ratio of upper reference limits
#' unless the cohort reports its own stratum.  Whether the measured baseline is raised is read
#' from the noise-free value.
#'
#' A patient whose baseline sits just below the limit reads positive at some test from the
#' assay's 9% imprecision alone, and over a year of tests almost surely, so equal-probability
#' slices alone place too little mass near the limit.  The slices are: n_lo of equal probability
#' in the bulk; thin slices at 4, 3, 2, 1 and 0.5 measurement-error widths below the limit and
#' 0.5, 1 and 2 above it, where the noise decides; and n_hi of equal probability in the upper
#' tail, where a background rate that rises with the baseline is decided.  Against a Monte Carlo
#' of single patients over s_b 0.5 to 2 and g_bg 0 and 2 the share positive is within 1.5%.
qsp2_pop_basenodes <- function(co, pp, n_lo = 4L, n_hi = 4L, thin_lo = c(4, 3, 2, 1, 0.5),
                               thin_hi = c(0.5, 1, 2)) {
  med <- pp$Tn0 * .assay_ratio(co); sdlog <- pp$s_b
  if (!is.null(co$base$sdlog)) sdlog <- co$base$sdlog
  if (!is.null(co$base$med)) med <- co$base$med
  if (!is.null(co$base$praised))
    med <- exp(log(co$limit) - sdlog * qnorm(co$base$praised, lower.tail = FALSE))
  zL <- (log(co$limit) - log(med)) / sdlog
  dz <- QSP2_SDLOG / sdlog                  # the assay's measurement error on this z scale
  a <- zL - max(thin_lo) * dz; h <- zL + max(thin_hi) * dz
  eq <- function(lo, hi, n) if (n < 2) numeric(0) else
    qnorm(pnorm(lo) + (pnorm(hi) - pnorm(lo)) * seq_len(n - 1) / n)
  cuts_lo <- c(-Inf, eq(-Inf, a, n_lo), zL - sort(thin_lo, decreasing = TRUE) * dz, zL)
  cuts_hi <- c(zL, zL + sort(thin_hi) * dz, eq(h, Inf, n_hi), Inf)
  slices <- function(cuts, raised) {
    lo <- head(cuts, -1); hi <- tail(cuts, -1)
    w <- pnorm(hi) - pnorm(lo); k <- w > 1e-300
    if (!any(k)) return(NULL)
    data.frame(z = ((dnorm(lo) - dnorm(hi)) / w)[k], w = w[k], raised = rep(raised, sum(k)))
  }
  out <- rbind(slices(cuts_lo, FALSE), slices(cuts_hi, TRUE))
  out$b <- med * exp(sdlog * out$z)
  if (co$rule == "normal_only") { out <- out[!out$raised, ]; out$w <- out$w / sum(out$w) }
  out
}

#' Probability that a scheduled test is positive from a background episode.
#' `L` is the threshold on the assay's own scale; `shape` is the time profile of the rate,
#' "constant" in the main analysis; "rising" or "falling" are alternatives, not used in the
#' calibration (both average 1 over the year).
#' @param b  the patient's own baseline troponin on the cohort's scale; the arrival rate is
#'   lam_bg * (b / Tn0)^g_bg, so a sicker patient has more background episodes.  Left NULL, or
#'   with g_bg zero or absent, the rate is the same for everyone.
qsp2_pop_pbg <- function(L, co, pp, day = NULL, shape = "constant", b = NULL) {
  # the background episodes are sized on the troponin T scale (C6); a troponin I threshold is
  # put on the reference troponin I scale by the cohort's assay ratio and then on
  # the troponin T scale by the ratio of the two reference limits
  L_T <- L * (if (co$kin == "TnT") 1 else QSP2_TNT_URL / (QSP2_REF_URL * .assay_ratio(co)))
  k_el <- if (co$kin == "TnT") QSP2_TROPONIN$hs_cTnT$k_el else QSP2_TROPONIN$hs_cTnI$k_el
  z <- (log(pp$m_bg) - log(L_T)) / pp$s_bg
  Dur <- ((log(pp$m_bg) - log(L_T)) * pnorm(z) + pp$s_bg * dnorm(z)) / k_el
  g <- if (is.null(day) || shape == "constant") 1 else if (shape == "rising")
    pmin(2, 2 * day / 365) else pmax(0, 2 * (1 - day / 365))
  gb <- if (is.null(pp$g_bg)) 0 else pp$g_bg
  rate <- pp$lam_bg * if (is.null(b) || gb == 0) 1 else
    (pmax(b, 1e-9) / (pp$Tn0 * .assay_ratio(co)))^gb
  1 - exp(-rate * g * Dur)
}

# --- the observation layer -------------------------------------------------------------------

QSP2_POP_SCALARS <- c("t2", "pos", "posdrug", "posbg", "dx", "dxs", "dxc", "dxpos", "dxpres",
                      "pres", "presraised", "sev", "sev0", "sevb40", "sevb30", "sevb0",
                      "o1num", "o1den", "o1bnum", "o3num", "o3den", "basep", "m1dx", "m2dx",
                      "ntests", "dmg", "dxcsev", "dxcnon", "posnopres")
QSP2_POP_MASSES  <- c("firstpos", "firstposdx", "dx", "dxc", "dxcsev", "dxcnon", "dmg", "lead")

.dayadd <- function(vec, days, v) {
  keep <- is.finite(days) & v != 0
  if (!any(keep)) return(vec)
  s <- rowsum(v[keep], days[keep], reorder = FALSE)
  i <- as.integer(rownames(s)) + 1L
  vec[i] <- vec[i] + s[, 1]
  vec
}

#' Everything a cohort observes, for every node of a trajectory grid.
#'
#' Each node is expanded into patient types, the presentation classes of qsp2_pop_prep times
#' 5 baseline points, all handled at once.  A test on the day of presentation is taken before
#' it.
#'
#' @param traj  from qsp2_pop_traj
#' @param prep  from qsp2_pop_prep, with the same observation parameters and follow-up
#' @param co    from qsp2_pop_cohort
#' @param no_treatment  TRUE for the untreated counterfactual (used by the checks)
#' @param sdlog  measurement error on the log scale; the assay's 9% by default, larger in a
#'               variant that adds within-person variation
#' @return list(S = node x scalar matrix, M = named list of node x day mass matrices,
#'              trop = node x bin matrix of troponin at clinical diagnosis, troppk = the same for
#'              the patient's highest troponin up to the halt of killing)
qsp2_pop_obs <- function(traj, prep, co, bg_shape = "constant", no_treatment = FALSE,
                         sdlog = QSP2_SDLOG) {
  pp <- prep$pp
  H <- min(traj$H, co$H); Hs <- traj$H
  if (prep$H != H) stop(sprintf("prep was made for %d days of follow-up, the cohort has %d", prep$H, H))
  N <- nrow(traj$nodes); n_day <- Hs + 1L
  sdl <- sdlog
  scale_a <- if (co$kin == "TnT") QSP2_TNT_RATIO else .assay_ratio(co)
  Ltx <- round(pp$L_tx)
  BN <- qsp2_pop_basenodes(co, pp); nB <- nrow(BN)
  Leff <- switch(co$rule,
    default     = ifelse(BN$raised, pmax(co$limit, co$r_rise * BN$b), co$limit),
    double      = pmax(co$abs_min, 2 * BN$b),
    conversion  = rep(co$limit, nB),
    normal_only = rep(co$limit, nB),
    threshold   = rep(co$limit, nB))
  can_pos <- as.numeric(!(co$rule == "conversion" & BN$raised))
  base0B  <- as.numeric(co$count_baseline & BN$raised)
  d <- co$days; nt <- length(d)
  pbgB <- if (nt) outer(seq_len(nB), d, function(i, dd)
    qsp2_pop_pbg(Leff[i], co, pp, dd, bg_shape, b = BN$b[i])) * can_pos else matrix(0, nB, 0)
  # testing exposure: pK[k + 1] = P(K = k) tests after the baseline, SK0[m + 1] =
  # P(K >= m), so test j is done with probability SK0[j + 1] unless a diagnosis came first
  pK <- qsp2_pop_expo_pk(co$expo, nt)
  SK0 <- rev(cumsum(rev(pK)))

  nM <- ncol(prep$TP); nR <- nM * nB
  mR <- rep(seq_len(nM), each = nB); bR <- rep(seq_len(nB), times = nM)
  forcedM <- prep$forced
  base0 <- base0B[bR]
  efmin <- QSP2_MYO$EFmin
  rr <- seq_len(nR)
  n_lt_day <- if (nt) findInterval((0:Hs) - 1, d) else rep(0L, Hs + 1L)   # tests before each day

  S  <- matrix(0, N, length(QSP2_POP_SCALARS), dimnames = list(NULL, QSP2_POP_SCALARS))
  M  <- lapply(setNames(QSP2_POP_MASSES, QSP2_POP_MASSES), function(z) matrix(0, N, n_day))
  TR <- TRP <- matrix(0, N, length(QSP2_TROP_BINS))
  one <- function(k) {
    wM <- prep$WM[k, ]; wR <- wM[mR] * BN$w[bR]
    loM <- prep$LO[k, ]; hiM <- prep$HI[k, ]
    Xa <- traj$X[k, ] * scale_a
    c2 <- traj$t2c[k]; t2 <- !is.na(c2) && c2 <= H
    gp <- prep$gpost[k, ]; cGH <- traj$cumG[k, H + 1L]
    tpM <- prep$TP[k, ]; tpM[!(wM > 0) | tpM > H] <- NA
    finM <- !is.na(tpM)
    clM <- if (t2) finM & (c2 <= tpM + QSP2_WORKUP) else rep(FALSE, nM)
    clM[is.na(clM)] <- FALSE
    tpR <- tpM[mR]; finR <- finM[mR]; clR <- clM[mR]
    tpRi <- ifelse(finR, tpR, Inf)
    # --- tests, and screening diagnosis at the first positive test after T2 and before
    #     a clinical diagnosis
    if (nt) {
      qdB <- pnorm((log(outer(BN$b, Xa[d + 1L], "+")) - log(Leff)) / sdl) * can_pos
      qB  <- 1 - (1 - qdB) * (1 - pbgB)
      q <- qB[bR, , drop = FALSE]
      # chronology: a patient diagnosed at presentation stops treatment, so the
      # tests scheduled after that day are not done (a test on the day comes before it)
      late <- outer(ifelse(clR, tpR, Inf), d, "<")
      q[late] <- 0
      SKx <- matrix(SK0[seq_len(nt) + 1L], nR, nt, byrow = TRUE)   # P(still on treatment at test j)
      tpe <- ifelse(clR, tpR, Inf)
      E <- if (t2) outer(tpe, d, ">") & matrix(d >= c2, nR, nt, byrow = TRUE) else
        matrix(FALSE, nR, nt)
      qE <- q * E
      CE <- matrixStats::rowCumprods(cbind(1, 1 - qE))
      # screening diagnosis at test j needs the test to be done; summing over the exposure,
      # P(no screening diagnosis) = 1 - sum_j P(K >= j) P(first eligible positive at j)
      Pj <- qE * CE[, seq_len(nt), drop = FALSE] * SKx; none_e <- 1 - rowSums(Pj)
      CQ <- matrixStats::rowCumprods(cbind(1, 1 - q))           # none positive among the first j
      # NQ[, m + 1]: no positive among the tests actually done up to test m, over the exposure:
      # sum_{k < m} P(K = k) CQ[, k + 1] + P(K >= m) CQ[, m + 1]
      NQ <- cbind(0, matrixStats::rowCumsums(CQ[, seq_len(nt), drop = FALSE] *
                                               matrix(pK[seq_len(nt)], nR, nt, byrow = TRUE))) +
        CQ * matrix(SK0, nR, nt + 1L, byrow = TRUE)
      # TL[, n + 1]: a positive among the tests done after test n (so K > n), by the recursion
      # U_n = P(K = n) + (1 - q_{n+1}) U_{n+1}, U_nt = P(K = nt); TL_n = P(K > n) - (1 - q_{n+1}) U_{n+1}
      U <- matrix(0, nR, nt + 1L); U[, nt + 1L] <- pK[nt + 1L]
      TL <- matrix(0, nR, nt + 1L)
      for (n in (nt - 1L):0) {
        U[, n + 1L] <- pK[n + 1L] + (1 - q[, n + 1L]) * U[, n + 2L]
        TL[, n + 1L] <- SK0[n + 2L] - (1 - q[, n + 1L]) * U[, n + 2L]
      }
    } else {
      q <- Pj <- SKx <- matrix(0, nR, 0); none_e <- rep(1, nR)
      CQ <- NQ <- CE <- matrix(1, nR, 1); TL <- matrix(0, nR, 1)
    }
    Pcl <- clR * none_e; Pscr <- rowSums(Pj); Pdx <- Pscr + Pcl; Pnone <- 1 - Pdx
    # --- severity in closed form over each class's critical values (T3 and the three bands),
    #     for the three courses: no diagnosis, clinical diagnosis, screening diagnosis at test j
    SN <- .sev_c(rep(cGH, nM), loM, hiM, forcedM)
    g_cl <- rep(1, nM)
    if (any(finM)) g_cl[finM] <- gp[pmin(tpM[finM] + Ltx, Hs) + 1L]
    if (no_treatment) g_cl[] <- cGH
    SC <- .sev_c(g_cl, loM, hiM, forcedM)
    if (nt) {
      g_t <- if (no_treatment) rep(cGH, nt) else gp[pmin(d + Ltx, Hs) + 1L]
      SS <- .sev_c(rep(g_t, each = nM), rep(loM, nt), rep(hiM, nt), FALSE)  # class varies fastest
      scr <- vapply(1:4, function(cc)
        rowSums(Pj * matrix(SS[, cc], nM, nt)[mR, , drop = FALSE]), numeric(nR))
    } else scr <- matrix(0, nR, 4)
    sevR <- if (t2) scr + Pcl * SC[mR, , drop = FALSE] + Pnone * SN[mR, , drop = FALSE] else
      matrix(0, nR, 4)
    sev0M <- if (t2) ifelse(matrix(clM, nM, 4), SC, SN) else matrix(0, nM, 4)
    # --- troponin measured at presentation
    pposR <- numeric(nR)
    if (any(finR)) {
      i <- which(finR)
      pposR[i] <- pnorm((log(BN$b[bR[i]] + Xa[tpR[i] + 1L]) - log(Leff[bR[i]])) / sdl) *
        can_pos[bR[i]]
    }
    # --- the first positive result: the baseline, then tests, then presentation on its day
    n_le_tp <- if (nt) findInterval(tpRi, d) else rep(0L, nR)
    if (nt) {
      A <- outer(tpRi, d, "<")                                  # a test after presentation
      fpT <- (1 - base0) * q * CQ[, seq_len(nt), drop = FALSE] * ifelse(A, 1 - pposR, 1) * SKx
    } else fpT <- matrix(0, nR, 0)
    fpP <- (1 - base0) * pposR * finR * NQ[cbind(rr, n_le_tp + 1L)]
    p_anypos <- 1 - (1 - base0) * NQ[, ncol(NQ)] * (1 - pposR * finR)
    # the chance of a diagnosis given that the first positive is that event (exact).  For a
    # first positive at presentation that does not lead to a clinical diagnosis, the diagnosis
    # needs a later positive test, and the patient must still be on treatment for it, so the
    # joint probability is the presentation term times CQ up to presentation times TL after it
    # (fpPdx); with full exposure this reduces to fpP * (1 - RQ).
    if (t2) {
      dxpB <- Pdx
      dxpT <- if (nt) ifelse(matrix(d < c2, nR, nt, byrow = TRUE), Pdx, 1) else fpT
      n_aft <- if (nt) findInterval(pmax(ifelse(finR, tpR, 0), c2 - 1), d) else rep(0L, nR)
      fpPdx <- ifelse(clR, fpP, (1 - base0) * pposR * finR * CQ[cbind(rr, n_le_tp + 1L)] *
                        TL[cbind(rr, n_aft + 1L)])
    } else { dxpB <- rep(0, nR); dxpT <- fpT * 0; fpPdx <- rep(0, nR) }
    # --- scalars
    s <- setNames(numeric(length(QSP2_POP_SCALARS)), QSP2_POP_SCALARS)
    s["t2"] <- as.numeric(t2)
    s["pos"] <- sum(wR * p_anypos)
    if (nt) {
      s["posdrug"] <- sum(wR * (1 - matrixStats::rowProds(1 - qdB[bR, , drop = FALSE])))
      s["posbg"]   <- sum(wR * (1 - matrixStats::rowProds(1 - pbgB[bR, , drop = FALSE])))
    }
    s["dx"] <- sum(wR * Pdx); s["dxs"] <- sum(wR * Pscr); s["dxc"] <- sum(wR * Pcl)
    n_lt_c2 <- if (nt && t2) findInterval(c2 - 1, d) else 0L
    s["dxpos"] <- sum(wR * (Pscr + Pcl * (1 - (1 - base0) * NQ[, n_lt_c2 + 1L] * (1 - pposR))))
    s["dxpres"] <- sum(wR * (Pcl + (if (nt) rowSums(Pj * outer(tpRi, d, "<=")) else 0)))
    s["pres"] <- sum(wR * finR)
    # a positive result in a patient who never presents within the follow-up; at a node that
    # reaches T2 this is a subclinical case found
    s["posnopres"] <- sum(wR * (1 - finR) * p_anypos)
    s["presraised"] <- sum(wR * Pcl * pposR)
    s["sev"] <- sum(wR * sevR[, 1])
    s[c("sevb40", "sevb30", "sevb0")] <- colSums(wR * sevR[, 2:4, drop = FALSE])
    s["sev0"] <- sum(wR * sev0M[mR, 1])
    s["basep"] <- sum(wR * BN$raised[bR])
    tcl <- ifelse(finR, tpR, 0)
    if (nt) {
      s["m1dx"] <- sum(wR * (as.vector(Pj %*% d) + Pcl * tcl))
      s["m2dx"] <- sum(wR * (as.vector(Pj %*% d^2) + Pcl * tcl^2))
      # a test is done if treatment still lasts (exposure), no screening diagnosis came at an
      # earlier test, and no clinical diagnosis came before its day (chronology)
      s["ntests"] <- sum(wR * rowSums(SKx * CE[, seq_len(nt), drop = FALSE] * !late))
    } else {
      s["m1dx"] <- sum(wR * Pcl * tcl); s["m2dx"] <- sum(wR * Pcl * tcl^2)
    }
    s["dxcsev"] <- sum(wR * Pcl * SC[mR, 1]) * t2
    s["dxcnon"] <- sum(wR * Pcl * (1 - SC[mR, 1])) * t2
    # O1: T2 cases that present within the follow-up; a positive result before presentation
    mlead <- numeric(n_day)
    if (t2 && any(finR)) {
      i <- which(finR)
      s["o1den"] <- sum(wR[i])
      n_lt_tp <- if (nt) findInterval(tpR[i] - 1, d) else rep(0L, length(i))
      s["o1num"] <- sum(wR[i] * (1 - (1 - base0[i]) * NQ[cbind(i, n_lt_tp + 1L)]))
      if (nt) {
        B <- outer(tpR[i], d, ">")
        s["o1bnum"] <- sum(wR[i] * rowSums(Pj[i, , drop = FALSE] * B))
        fb <- (1 - base0[i]) * q[i, , drop = FALSE] * CQ[i, seq_len(nt), drop = FALSE] * B *
          SKx[i, , drop = FALSE]
        lag <- outer(tpR[i], d, "-")
        mlead <- .dayadd(mlead, lag[B], (wR[i] * fb)[B])
      }
      mlead <- .dayadd(mlead, tpR[i], wR[i] * base0[i])
    }
    # O3: patients severe under clinical diagnosis alone; a positive result before their LVEF
    # meets the criterion.  Presenters by LVEF cross on their presentation day.  For the
    # others, by day t the patients of a class who have crossed are those whose critical value
    # exceeds cumG(t); in a clinically diagnosed course the rest of the T3 mass is placed on the
    # halt, after which the course is no longer the untreated one.
    sev0R <- sev0M[mR, 1]
    if (t2 && any(sev0M[, 1] > 0)) {
      s["o3den"] <- sum(wR * sev0R)
      cg <- traj$cumG[k, seq_len(H + 1L)]; sv <- 0:H
      FM <- matrix(0, nM, H + 1L)
      for (m in which(sev0M[, 1] > 0 & wM > 0)) {
        if (forcedM[m]) { if (finM[m]) FM[m, ] <- as.numeric(sv >= tpM[m]); next }
        pset <- qsp2_S_c(loM[m]) - qsp2_S_c(hiM[m])
        if (!(pset > 0)) next
        if (clM[m] && !no_treatment) {
          th0 <- min(tpM[m] + Ltx, H)
          f <- pmax(qsp2_S_c(pmax(cg, g_cl[m], loM[m])) - qsp2_S_c(hiM[m]), 0) / pset
          f[sv >= th0] <- sev0M[m, 1]
        } else f <- pmax(qsp2_S_c(pmax(cg, loM[m])) - qsp2_S_c(hiM[m]), 0) / pset
        FM[m, ] <- f
      }
      dF <- FM - cbind(0, FM[, -ncol(FM), drop = FALSE])
      Pb <- 1 - (1 - base0) * NQ[, n_lt_day[seq_len(H + 1L)] + 1L, drop = FALSE]
      s["o3num"] <- sum(wR * rowSums(dF[mR, , drop = FALSE] * Pb))
    }
    # --- masses over days
    mfp <- mfpd <- mdx <- mdxc <- mdxs <- mdxn <- mdmg <- numeric(n_day)
    mfp  <- .dayadd(mfp, rep(0, nR), wR * base0)
    mfpd <- .dayadd(mfpd, rep(0, nR), wR * base0 * dxpB)
    if (nt) {
      dd <- rep(d, each = nR)
      mfp  <- .dayadd(mfp, dd, as.vector(wR * fpT))
      mfpd <- .dayadd(mfpd, dd, as.vector(wR * fpT * dxpT))
      mdx  <- .dayadd(mdx, dd, as.vector(wR * Pj))
      # C3 damage: a positive result of any cause, or a fall in LVEF of more than 10 points to
      # below 53% at a scheduled test (Nishikawa 2022's definition).  The LVEF arm
      # holds for baselines in (EFmin + 10 / (1 - g), EFmin + 43 / g), with g the lowest G seen
      # at a test so far, and it is integrated over each class's baselines.
      gmin <- cummin(traj$G[k, d + 1L])
      JL <- ifelse(gmin < 1, efmin + QSP2_EF_FALL / (1 - gmin), Inf)
      JU <- efmin + (QSP2_C3_EF - efmin) / pmax(gmin, 1e-300)
      pset <- qsp2_S_c(loM) - qsp2_S_c(hiM)
      Ff <- matrix(.p_in_set(rep(JL, each = nM), rep(JU, each = nM), rep(loM, nt), rep(hiM, nt)),
                   nM, nt) / ifelse(pset > 0, pset, 1)
      Ff <- pmin(Ff, 1)[mR, , drop = FALSE]
      Fprev <- cbind(0, Ff[, -nt, drop = FALSE])
      # a scheduled visit is only there while treatment lasts and before a clinical diagnosis
      dm <- CQ[, seq_len(nt), drop = FALSE] * pmax((1 - Fprev) - (1 - Ff) * (1 - q), 0) * SKx * !late
      mdmg <- .dayadd(mdmg, dd, as.vector(wR * dm))
      s["dmg"] <- sum(wR * rowSums(dm))
    }
    mfp  <- .dayadd(mfp, tpR, wR * fpP)
    mfpd <- .dayadd(mfpd, tpR, wR * fpPdx)
    tr <- trp <- numeric(length(QSP2_TROP_BINS))
    if (any(clR)) {
      i <- which(clR)
      mdx  <- .dayadd(mdx, tpR[i], wR[i] * Pcl[i])
      mdxc <- .dayadd(mdxc, tpR[i], wR[i] * Pcl[i])
      mdxs <- .dayadd(mdxs, tpR[i], wR[i] * Pcl[i] * SC[mR[i], 1])
      mdxn <- .dayadd(mdxn, tpR[i], wR[i] * Pcl[i] * (1 - SC[mR[i], 1]))
      lv <- log10(pmax(BN$b[bR[i]] + Xa[tpR[i] + 1L], 1))
      tr <- .dayadd(tr, findInterval(lv, QSP2_TROP_BINS, all.inside = TRUE) - 1L, wR[i] * Pcl[i])
      # the patient's highest troponin, which is what a source reporting the highest value during
      # therapy measures.  Killing halts L_tx days after diagnosis and troponin falls after that,
      # so the highest value is reached by then.
      cmx <- cummax(Xa)
      pk <- log10(pmax(BN$b[bR[i]] + cmx[pmin(tpR[i] + Ltx, Hs) + 1L], 1))
      trp <- .dayadd(trp, findInterval(pk, QSP2_TROP_BINS, all.inside = TRUE) - 1L, wR[i] * Pcl[i])
    }
    list(s = s, m = list(firstpos = mfp, firstposdx = mfpd, dx = mdx, dxc = mdxc,
                         dxcsev = mdxs, dxcnon = mdxn, dmg = mdmg, lead = mlead),
         tr = tr, trp = trp)
  }
  # Nodes where nothing happens (no drug-released troponin to speak of, no presentation, no
  # T2, no loss of function) all observe the same thing, so that is computed once.
  pk <- if (is.null(traj$peakX)) apply(traj$X, 1, max) else traj$peakX
  pres_w <- rowSums(prep$WM[, prep$cls != "N", drop = FALSE])
  null <- traj$ok & pk < 1e-4 & is.na(traj$t2c) & pres_w == 0 &
    apply(traj$G, 1, function(g) all(g > 1 - 1e-9))
  null[is.na(null)] <- FALSE
  done <- NULL
  for (k in which(traj$ok)) {
    r <- if (null[k] && !is.null(done)) done else one(k)
    if (null[k] && is.null(done)) done <- r
    S[k, ] <- r$s
    for (nm in QSP2_POP_MASSES) M[[nm]][k, ] <- r$m[[nm]]
    TR[k, ] <- r$tr; TRP[k, ] <- r$trp
  }
  list(S = S, M = M, trop = TR, troppk = TRP, cohort = co, H = H, n_null = sum(null))
}

#' Observe a cohort that may be a mixture of calendars (qsp2_pop_cohort_mix): the weighted sum of
#' its parts' observations, every one of which is linear in the patients.
qsp2_pop_obs_any <- function(traj, prep, co, ...) {
  if (is.null(co$parts)) return(qsp2_pop_obs(traj, prep, co, ...))
  os <- lapply(co$parts, function(p) qsp2_pop_obs(traj, prep, p, ...))
  w <- co$w_parts
  add <- function(get) Reduce(`+`, Map(function(o, wi) wi * get(o), os, w))
  list(S = add(function(o) o$S),
       M = setNames(lapply(names(os[[1]]$M), function(nm) add(function(o) o$M[[nm]])), names(os[[1]]$M)),
       trop = add(function(o) o$trop), troppk = add(function(o) o$troppk),
       cohort = co, H = os[[1]]$H, n_null = os[[1]]$n_null)
}

#' Mix the two regimens' observations for a cohort's own share of combination therapy.
qsp2_pop_mix <- function(o_mono, o_combo, w) {
  if (w <= 0) return(o_mono)
  if (w >= 1) return(o_combo)
  list(S = (1 - w) * o_mono$S + w * o_combo$S,
       M = Map(function(a, b) (1 - w) * a + w * b, o_mono$M, o_combo$M),
       trop = (1 - w) * o_mono$trop + w * o_combo$trop,
       troppk = (1 - w) * o_mono$troppk + w * o_combo$troppk,
       cohort = o_mono$cohort, H = o_mono$H)
}

# --- population weights ---------------------------------------------------------------------

#' Probability mass of each grid cell for a lognormal, cells bounded at geometric midpoints.
qsp2_pop_cellmass <- function(grid, meanlog, sdlog) {
  lg <- log(grid); e <- c(-Inf, (head(lg, -1) + tail(lg, -1)) / 2, Inf)
  diff(pnorm(e, mean = meanlog, sd = sdlog))
}

#' Node weights for many population triples at once: a (triples x nodes) matrix.
#' Mass on nodes whose solve failed is dropped and the rest renormalised; the dropped mass is
#' returned as attribute "lost" so a caller can refuse a set that loses too much.
qsp2_pop_W <- function(traj, tri) {
  nodes <- traj$nodes
  WS <- vapply(seq_len(nrow(tri)), function(i)
    qsp2_pop_cellmass(traj$sgrid, tri$mu[i], tri$sigma[i]), numeric(length(traj$sgrid)))
  WT <- vapply(seq_len(nrow(tri)), function(i)
    qsp2_pop_cellmass(traj$taugrid, 0, tri$sd_tau[i]), numeric(length(traj$taugrid)))
  W <- t(matrix(WS[nodes$i_s, ], ncol = nrow(tri)) * matrix(WT[nodes$i_t, ], ncol = nrow(tri)))
  lost <- rowSums(W[, !traj$ok, drop = FALSE])
  W[, !traj$ok] <- 0
  W <- W / pmax(rowSums(W), 1e-300)
  attr(W, "lost") <- lost
  W
}

#' Quantiles of a population time distribution from a cumulative mass matrix (rows are
#' populations, columns days 0..H).  NA where the event has no mass.  The mass on day t is
#' spread evenly over (t - 0.5, t + 0.5], so an event on day 56 alone has median 56.
qsp2_pop_tq_cum <- function(C, probs) {
  C <- as.matrix(C)
  tot <- C[, ncol(C)]
  out <- matrix(NA_real_, nrow(C), length(probs))
  rows <- seq_len(nrow(C))
  for (j in seq_along(probs)) {
    tgt <- probs[j] * tot
    idx <- max.col((C >= tgt - 1e-12 * tot) * 1, ties.method = "first")
    prev <- ifelse(idx > 1, C[cbind(rows, pmax(idx - 1, 1))], 0)
    cur  <- C[cbind(rows, idx)]
    frac <- ifelse(cur > prev, (tgt - prev) / (cur - prev), 1)
    out[, j] <- ifelse(tot > 0, pmax(idx - 1.5 + frac, 0), NA_real_)
  }
  out
}

#' Quantiles of a population time distribution, from node masses, for every triple.
qsp2_pop_tq <- function(W, mass, probs)
  qsp2_pop_tq_cum(W %*% matrixStats::rowCumsums(mass), probs)

#' Weighted quantile of a node-level value (e.g. t1ref) among nodes, for every triple.
qsp2_pop_nodeq <- function(W, value, wnode, probs) {
  ok <- is.finite(value) & wnode > 0
  if (!any(ok)) return(matrix(NA_real_, nrow(W), length(probs)))
  o <- order(value[ok]); v <- value[ok][o]
  Wn <- sweep(W[, ok, drop = FALSE][, o, drop = FALSE], 2, wnode[ok][o], "*")
  C <- matrixStats::rowCumsums(Wn); tot <- C[, ncol(C)]
  out <- sapply(probs, function(pq) {
    idx <- max.col((C >= pq * tot - 1e-12 * tot) * 1, ties.method = "first")
    ifelse(tot > 0, v[idx], NA_real_)
  })
  matrix(out, nrow(W))
}

#' Median troponin (ng/L) from population histograms over QSP2_TROP_BINS (rows are
#' populations), interpolated within the 0.05-log10 bins.
qsp2_pop_tropmed_hist <- function(Hm) {
  C <- matrixStats::rowCumsums(as.matrix(Hm)); tot <- C[, ncol(C)]
  rows <- seq_len(nrow(C))
  idx <- max.col((C >= 0.5 * tot) * 1, ties.method = "first")
  prev <- ifelse(idx > 1, C[cbind(rows, pmax(idx - 1, 1))], 0)
  cur  <- C[cbind(rows, idx)]
  frac <- ifelse(cur > prev, (0.5 * tot - prev) / (cur - prev), 0.5)
  ifelse(tot > 0, 10^(QSP2_TROP_BINS[idx] + 0.05 * frac), NA_real_)
}

#' Median of the population's troponin at clinical diagnosis (ng/L), for every triple.
qsp2_pop_tropmed <- function(W, trop) qsp2_pop_tropmed_hist(W %*% trop)

#' The observation parameters of a parameter vector, as the engine reads them.
qsp2_pop_pp <- function(pars)
  as.list(pars[c("k_sym", "s_ksym", "Tn0", "s_b", "L_tx", "lam_bg", "m_bg", "s_bg",
                 "r_rep", "r_scar", "g_bg")])
