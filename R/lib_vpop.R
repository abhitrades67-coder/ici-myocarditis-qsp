# ---------------------------------------------------------------------------
# lib_vpop.R  -  virtual populations, built from the susceptibility response curve.
#
# Every outcome the model produces is monotone in the patient's precursor supply S:
# a larger S gives a larger infiltrate, more damage and an earlier crossing of each
# tier.  So a whole virtual population can be summarised by running the model once
# over a grid of S and then integrating over the distribution of S, instead of
# simulating thousands of patients one at a time.  That is what makes calibration
# affordable; the monotonicity it relies on is checked, not assumed.
#
# Definitions only.
# ---------------------------------------------------------------------------

QSP2_SGRID <- 10^seq(-4, 4, by = 0.2)     # 41 points, five orders of magnitude

#' Run one parameter set over the susceptibility grid and record, for each S, whether
#' each tier is reached and when.
#'
#' @return a data.table with one row per S
qsp2_response_curve <- function(pars, regimen, times = seq(0, 365, by = 1),
                                assay = "hs_cTnI", sgrid = QSP2_SGRID,
                                schedules = c("S0", "S1", "S2", "S3", "S4"),
                                extra_schedules = NULL, cohorts = NULL) {
  lim <- QSP2_TROPONIN[[assay]]$limit
  t2_cells <- qsp2_t2_cells(pars)
  sch <- lapply(schedules, function(s) qsp2_schedule(s, regimen, max(times)))
  names(sch) <- schedules
  # a cohort's own testing calendar, when it is not one of the five compared schedules.
  # These stay on the reference assay and threshold; a cohort passed through `cohorts`
  # carries its own.
  if (length(extra_schedules)) sch <- c(sch, extra_schedules)
  # Simulate once per assay that is actually asked for, not once per cohort.
  alt <- setdiff(unique(vapply(cohorts, function(x) x$assay, character(1))), assay)
  rows <- lapply(sgrid, function(s) {
    p <- pars; p[["Sus"]] <- s
    r <- tryCatch(qsp2_simulate(p, regimen, times, assay = assay), error = function(e) NULL)
    if (is.null(r) || anyNA(r$Tn) || anyNA(r$v)) return(NULL)
    ra <- lapply(alt, function(a)
      tryCatch(qsp2_simulate(p, regimen, times, assay = a), error = function(e) NULL))
    names(ra) <- alt
    if (length(ra) && any(vapply(ra, is.null, logical(1)))) return(NULL)
    ra[[assay]] <- r
    ef <- qsp2_lvef(r$v, r$d, p)
    # tier T1: troponin above the assay limit
    w1 <- which(r$Tn > lim)
    # tier T2: infiltrate above the myocarditis threshold, with injury
    w2 <- which(r$Tcells_heart > t2_cells & r$Tn > lim)
    # tier T3: ejection fraction below the registry cut-off
    w3 <- which(ef < 50)
    # Clinical presentation: a troponin high enough to be acted on, or an
    # ejection fraction already below the cut-off.  Presentation in this disease follows
    # acute inflammation, which troponin tracks, not the slow accumulation of lost tissue.
    wc <- which(r$Tn > p[["k_sym"]] * lim | ef < 50)
    out <- data.table(S = s,
               T1 = length(w1) > 0, t1 = if (length(w1)) r$t[w1[1]] else NA_real_,
               T2 = length(w2) > 0, t2 = if (length(w2)) r$t[w2[1]] else NA_real_,
               T3 = length(w3) > 0, t3 = if (length(w3)) r$t[w3[1]] else NA_real_,
               CL = length(wc) > 0, tc = if (length(wc)) r$t[wc[1]] else NA_real_,
               peakTn = max(r$Tn), loss = 1 - min(r$v), lvef = min(ef),
               cells = max(r$Tcells_heart),
               # the troponin at the moment of presentation, which is not always above
               # the limit: presentation can also be triggered by the ejection fraction
               clinTn = if (length(wc)) r$Tn[wc[1]] else NA_real_)
    # when each schedule would first see a raised result, and how many tests it costs
    for (nm in names(sch)) {
      d <- sch[[nm]]
      j <- match(d, r$t)
      j <- j[!is.na(j)]
      hit <- j[r$Tn[j] > lim]
      out[[paste0("det_", nm)]] <- length(hit) > 0
      out[[paste0("tdet_", nm)]] <- if (length(hit)) r$t[hit[1]] else NA_real_
      out[[paste0("ntest_", nm)]] <- length(j)
    }
    # the same, for a cohort scored at its own assay and its own threshold
    for (nm in names(cohorts)) {
      co <- cohorts[[nm]]
      rc <- ra[[co$assay]]
      j <- match(co$days, rc$t)
      j <- j[!is.na(j)]
      hit <- j[rc$Tn[j] > co$limit]
      out[[paste0("det_", nm)]] <- length(hit) > 0
      out[[paste0("tdet_", nm)]] <- if (length(hit)) rc$t[hit[1]] else NA_real_
      out[[paste0("ntest_", nm)]] <- length(j)
    }
    out
  })
  ok <- !vapply(rows, is.null, logical(1))
  if (!any(ok)) return(NULL)
  rbindlist(rows[ok])
}

#' Is the curve monotone in S, as the virtual-population shortcut assumes?
#' Returns the share of adjacent pairs that respect monotonicity for each outcome.
qsp2_check_monotone <- function(cur) {
  mono <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(1); mean(diff(x) >= -1e-9) }
  anti <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(1); mean(diff(x) <= 1e-9) }
  c(loss = mono(cur$loss), peakTn = mono(cur$peakTn), cells = mono(cur$cells),
    lvef = anti(cur$lvef), t1 = anti(cur$t1), t2 = anti(cur$t2), tc = anti(cur$tc))
}

#' The susceptibility at which a tier is first reached, read from the curve by
#' interpolating on the log scale between the two bracketing grid points.
qsp2_threshold <- function(cur, tier) {
  h <- cur[[tier]]
  k <- which(h)
  if (!length(k)) return(Inf)
  if (k[1] == 1L) return(cur$S[1])
  # linear interpolation in log S is not possible for a step, so take the geometric
  # midpoint of the bracketing points; the grid is fine enough that this is small
  sqrt(cur$S[k[1] - 1L] * cur$S[k[1]])
}

#' Share of a lognormal population above a threshold.
qsp2_share_above <- function(thr, mu, sigma) {
  if (!is.finite(thr)) return(0)
  stats::pnorm(log(thr), mean = mu, sd = sigma, lower.tail = FALSE)
}

#' Quantiles of the time at which a tier is crossed, over the patients who cross it.
#' The time is read from the curve and weighted by the lognormal density of S.
qsp2_time_quantiles <- function(cur, tcol, hcol, mu, sigma, probs = c(0.25, 0.5, 0.75)) {
  ok <- cur[[hcol]] & is.finite(cur[[tcol]])
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  s <- cur$S[ok]; tt <- cur[[tcol]][ok]
  # weight each grid point by the probability mass of S it represents
  ls <- log(s); edges <- c(-Inf, (head(ls, -1) + tail(ls, -1)) / 2, Inf)
  w <- diff(stats::pnorm(edges, mean = mu, sd = sigma))
  if (sum(w) <= 0) return(rep(NA_real_, length(probs)))
  o <- order(tt); tt <- tt[o]; w <- w[o] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) tt[which(cw >= p)[1]], numeric(1))
}

#' Everything a calibration target needs, for one parameter set and one regimen.
qsp2_population_summary <- function(cur, mu, sigma) {
  list(
    p_T1   = qsp2_share_above(qsp2_threshold(cur, "T1"), mu, sigma),
    p_T2   = qsp2_share_above(qsp2_threshold(cur, "T2"), mu, sigma),
    p_T3   = qsp2_share_above(qsp2_threshold(cur, "T3"), mu, sigma),
    p_clin = qsp2_share_above(qsp2_threshold(cur, "CL"), mu, sigma),
    t_clin = qsp2_time_quantiles(cur, "tc", "CL", mu, sigma),
    t_T1   = qsp2_time_quantiles(cur, "t1", "T1", mu, sigma),
    t_T2   = qsp2_time_quantiles(cur, "t2", "T2", mu, sigma))
}


# --- background troponin rises (Section 8 of the specification) -------------------------
# Rises from causes other than this drug reaction are what make the frequency of testing
# matter: a schedule that tests more often catches more of them, and they are the source
# of false alarms.  They are modelled where the specification puts them, in the
# observation model rather than in the differential equations: each scheduled test has a
# probability p_bg of returning a raised result for some other reason, independently of
# the drug.  A schedule of n tests therefore returns at least one such result with
# probability 1 - (1 - p_bg)^n.
#
# p_bg is calibrated to C6: in 164 patients screened prospectively, 26 had a significant
# rise and 8 of those were myocarditis, so 18 of 164 had a rise from another cause.

#' Probability that a schedule of n tests sees at least one background rise.
qsp2_p_background <- function(p_bg, n_tests) 1 - (1 - p_bg)^pmax(n_tests, 0)

#' Probability that a schedule returns a raised result from ANY cause: the drug-driven
#' rise or a background rise.
#'
#' Any block that needs this quantity calls this function, so there is only one definition
#' of it.
qsp2_p_raised <- function(p_drug, p_bg, n_tests)
  1 - (1 - p_drug) * (1 - qsp2_p_background(p_bg, n_tests))

#' Detection under one schedule, combining the drug-driven rise with background rises.
#'
#' @return a list with the share detected, the share detected who have myocarditis, the
#'   false alarms per case found, and the positive predictive value.
qsp2_detection <- function(cur, sched, mu, sigma, p_bg) {
  n   <- cur[[paste0("ntest_", sched)]][1]
  p_d <- qsp2_share_above(qsp2_threshold(cur, paste0("det_", sched)), mu, sigma)
  p_any <- qsp2_p_raised(p_d, p_bg, n)
  # a patient who reaches the myocarditis tier is above the detection threshold, so the
  # schedule finds them whenever it finds anyone driven by the drug
  p_T2  <- qsp2_share_above(qsp2_threshold(cur, "T2"), mu, sigma)
  p_T2_found <- min(p_T2, p_d)
  list(p_detected = p_any, p_found_T2 = p_T2_found,
       false_per_case = if (p_T2_found > 0) max(0, p_any - p_T2_found) / p_T2_found else NA_real_,
       ppv = if (p_any > 0) p_T2_found / p_any else NA_real_,
       n_tests = n)
}
