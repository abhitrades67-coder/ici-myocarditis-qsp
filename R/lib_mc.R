# ---------------------------------------------------------------------------
# lib_mc.R  -  a patient-by-patient Monte Carlo, written separately from lib_pop.R.
#
# It exists to check the population engine.  It draws each patient's grid node, symptom
# threshold, baseline LVEF, baseline troponin, the measurement error of every test and the
# background episodes (as an actual Poisson process of decaying peaks), and applies the rules
# of the specification directly to each patient's LVEF and troponin.
# Where the engine uses closed forms (severity within LVEF slices, the post-halt nadir, the
# background probability), this file does not: it follows the treated course explicitly and
# simulates the episodes.  The only things shared with the engine are the stored trajectories
# and the conventions of the specification.
#
# Definitions only: sourcing this file runs nothing.
# ---------------------------------------------------------------------------

#' Simulate n patients of one population (one row of `tri`) through one cohort.
#' With `kk` given, patient i sits on node kk[i] and n is length(kk), so that each patient can
#' sit on a node solved with its own parameters.
#' Returns one row per patient.
qsp2_mc_cohort <- function(traj, pp, co, tri, n, seed = 1, kk = NULL) {
  set.seed(seed)
  pp <- as.list(pp)
  pars <- traj$pars
  if (is.null(kk)) {
    W <- qsp2_pop_W(traj, tri)[1, ]
    kk <- sample.int(length(W), n, replace = TRUE, prob = W)
  } else n <- length(kk)
  H <- min(traj$H, co$H); Hs <- traj$H
  em <- QSP2_MYO$EFmin
  # a troponin I assay other than the reference sits on the
  # reference scale times the ratio of the limits, or (conv "none") on the reference scale itself
  r_assay <- if (co$kin == "TnT") co$url / QSP2_REF_URL else
    if (identical(co$conv, "none")) 1 else co$url / QSP2_REF_URL
  scale_a <- if (co$kin == "TnT") QSP2_TNT_RATIO else r_assay
  kel <- if (co$kin == "TnT") QSP2_TROPONIN$hs_cTnT$k_el else QSP2_TROPONIN$hs_cTnI$k_el
  ksym <- pp$k_sym * exp(pp$s_ksym * rnorm(n))
  ef0 <- QSP2_MYO$EF0 + QSP2_EF_SD * rnorm(n)
  med <- pp$Tn0 * r_assay; sdb <- pp$s_b
  if (!is.null(co$base$sdlog)) sdb <- co$base$sdlog
  if (!is.null(co$base$med)) med <- co$base$med
  if (!is.null(co$base$praised))
    med <- exp(log(co$limit) - sdb * qnorm(co$base$praised, lower.tail = FALSE))
  b <- med * exp(sdb * rnorm(n)); raised <- b > co$limit
  keep <- if (co$rule == "normal_only") !raised else rep(TRUE, n)
  Leff <- switch(co$rule,
    default = ifelse(raised, pmax(co$limit, co$r_rise * b), co$limit),
    double = pmax(co$abs_min, 2 * b), rep(co$limit, n))
  canp <- !(co$rule == "conversion" & raised)
  d <- co$days; nt <- length(d)
  sdl <- sqrt(log(1 + QSP2_CV_ANA^2))
  Ltx <- round(pp$L_tx)
  rr <- pp$r_rep + pp$r_scar
  # each patient's number of tests after the baseline: drawn once, before the
  # course is known; every scheduled test when the cohort has no exposure model
  pK <- qsp2_pop_expo_pk(co$expo, nt)
  Kx <- sample.int(nt + 1L, n, replace = TRUE, prob = pK) - 1L
  # background: Poisson episodes over the 6 days before a test, lognormal peaks on the
  # troponin T scale, decaying at the assay's clearance rate.  The arrival rate is
  # lam_bg (b / reference baseline)^g_bg for a patient with baseline b.
  gb <- if (is.null(pp$g_bg)) 0 else pp$g_bg
  bg <- function(L, bi) {
    L_T <- L * (if (co$kin == "TnT") 1 else QSP2_TNT_URL / (QSP2_REF_URL * r_assay))
    rate <- pp$lam_bg * (if (gb == 0) 1 else (pmax(bi, 1e-9) / (pp$Tn0 * r_assay))^gb)
    ne <- rpois(length(L), rate * 6)
    hit <- logical(length(L))
    for (i in which(ne > 0)) {
      age <- runif(ne[i], 0, 6); M <- pp$m_bg * exp(pp$s_bg * rnorm(ne[i]))
      hit[i] <- any(M * exp(-kel * age) > L_T[i])
    }
    hit
  }
  # LVEF criterion for baseline e at G: below 50% and at least 10 points below baseline
  crit <- function(e, g) { ef <- em + (e - em) * g; ef < QSP2_EF_CUT & e - ef >= QSP2_EF_FALL }
  # the nadir of G within the follow-up for a course whose killing halts on day th, following
  # the treated course explicitly (repair and scarring after the halt) up to H
  nadir_treated <- function(k, th) {
    if (th >= H) return(min(traj$G[k, seq_len(H + 1L)]))
    g_before <- min(traj$G[k, seq_len(th + 1L)])
    vh <- traj$V[k, th + 1L]; dh <- traj$D[k, th + 1L]
    s <- seq(0, H - th, by = 0.25)
    v <- vh + dh * pp$r_rep / rr * (1 - exp(-rr * s)); dd <- dh * exp(-rr * s)
    # treatment resolves the infiltrate as well as stopping killing, so the
    # post-halt course is followed on tissue alone
    min(g_before, min((v + pars[["eta"]] * dd)^pars[["pexp"]]))
  }
  out <- data.frame(keep = keep, raised = raised, t2 = FALSE, pos = FALSE, dx = FALSE, dxs = FALSE,
                    dxc = FALSE, sev = FALSE, sev0 = FALSE, o1d = FALSE, o1n = FALSE, o3d = FALSE,
                    o3n = FALSE, dmg = FALSE, dxpos = FALSE, dxpres = FALSE, presr = FALSE,
                    tdx = NA_real_, tfp = NA_real_, tdmg = NA_real_, tp = NA_real_)
  for (k in unique(kk)) {
    ii <- which(kk == k); m <- length(ii)
    X <- traj$X[k, ]; Xa <- X * scale_a; G <- traj$G[k, ]
    c2 <- traj$t2c[k]; t2 <- !is.na(c2) && c2 <= H
    P <- matrix(FALSE, m, max(nt, 1))
    if (nt) for (j in seq_len(nt)) {
      meas <- (b[ii] + Xa[d[j] + 1]) * exp(sdl * rnorm(m))
      P[, j] <- canp[ii] & (meas > Leff[ii] | bg(Leff[ii], b[ii]))
    }
    for (a in seq_len(m)) {
      i <- ii[a]
      trig <- which((X + pp$Tn0 > max(ksym[i], 1) * QSP2_REF_URL) | crit(ef0[i], G))
      tp <- if (length(trig) && trig[1] - 1 <= H) trig[1] - 1 else NA
      cl <- t2 && !is.na(tp) && c2 <= tp + QSP2_WORKUP
      ep <- rnorm(1)          # drawn for every patient, so paired runs keep the same random stream
      ppos <- !is.na(tp) && canp[i] && (b[i] + Xa[tp + 1]) * exp(sdl * ep) > Leff[i]
      # a scheduled test is done only while treatment lasts (the patient's K) and not after a
      # clinical diagnosis, which stops treatment
      perf <- if (nt) seq_len(nt) <= Kx[i] & (if (cl) d <= tp else rep(TRUE, nt)) else logical(0)
      pj <- if (nt) P[a, seq_len(nt)] & perf else logical(0)
      tpe <- if (cl) tp else Inf
      scr <- NA
      if (t2 && nt) { e <- which(d >= c2 & d < tpe & pj); if (length(e)) scr <- d[e[1]] }
      tdx <- if (!is.na(scr)) scr else if (cl) tp else NA
      g_star <- if (is.na(tdx)) min(G[seq_len(H + 1L)]) else nadir_treated(k, tdx + Ltx)
      sev <- t2 && crit(ef0[i], g_star)
      g0 <- if (cl) nadir_treated(k, tp + Ltx) else min(G[seq_len(H + 1L)])
      sev0 <- t2 && crit(ef0[i], g0)
      base0 <- co$count_baseline && raised[i]
      ev <- c(if (base0) 0, if (nt) d[pj], if (!is.na(tp) && ppos) tp + 0.5)
      o1d <- t2 && !is.na(tp)
      # O3: the day the criterion is met on the S0 course, the untreated path up to the halt
      # and the halt itself for later crossings
      th0 <- if (cl) min(tp + Ltx, H) else H + 1
      z <- which(crit(ef0[i], G[seq_len(H + 1L)]) & (0:H) < th0)
      t50 <- if (length(z)) z[1] - 1 else th0
      evt <- c(if (base0) 0, if (nt) d[pj])
      # C3 damage at a scheduled test: a positive result, or a fall of more than 10 points to
      # below 53%
      if (nt) {
        eft <- em + (ef0[i] - em) * G[d + 1]
        f <- perf & (pj | (eft < QSP2_C3_EF & ef0[i] - eft > QSP2_EF_FALL))
        if (any(f)) { out$dmg[i] <- TRUE; out$tdmg[i] <- d[which(f)[1]] }
      }
      dxpos <- !is.na(scr) || (cl && is.na(scr) &&
        (base0 || (nt && any(pj & d < c2)) || ppos))
      out[i, c("t2", "pos", "dx", "dxs", "dxc", "sev", "sev0", "o1d", "o1n", "o3d", "o3n",
               "dxpos", "dxpres", "presr")] <-
        list(t2, length(ev) > 0, !is.na(tdx), !is.na(scr), is.na(scr) && cl, sev, sev0,
             o1d, o1d && length(ev) > 0 && any(ev < tp), sev0,
             sev0 && length(evt) > 0 && any(evt < t50), !is.na(tdx) && dxpos,
             !is.na(tdx) && !is.na(tp) && tp <= tdx, is.na(scr) && cl && ppos)
      out$tdx[i] <- tdx
      out$tfp[i] <- if (length(ev)) min(ev) else NA
      out$tp[i] <- tp
    }
  }
  out[out$keep, ]
}

#' The statistics the engine and the Monte Carlo are compared on, with the Monte Carlo's own
#' standard errors.
qsp2_mc_summary <- function(sm) {
  sh <- function(v) c(mean(v), sqrt(mean(v) * (1 - mean(v)) / length(v)))
  ratio <- function(num, den) { p <- sum(num) / sum(den); c(p, sqrt(p * (1 - p) / max(sum(den), 1))) }
  rbind(t2 = sh(sm$t2), pos = sh(sm$pos), dx = sh(sm$dx), dxs = sh(sm$dxs), dxc = sh(sm$dxc),
        sev = sh(sm$sev), sev0 = sh(sm$sev0), O1 = ratio(sm$o1n, sm$o1d),
        O3 = ratio(sm$o3n, sm$o3d), basep = sh(sm$raised), dmg = sh(sm$dmg),
        dxpos = sh(sm$dxpos), dxpres = sh(sm$dxpres), presraised = sh(sm$presr),
        mean_tdx = c(mean(sm$tdx, na.rm = TRUE), sd(sm$tdx, na.rm = TRUE) / sqrt(sum(!is.na(sm$tdx)))),
        med_tdx = c(median(sm$tdx, na.rm = TRUE), NA),
        med_tfp = c(median(sm$tfp, na.rm = TRUE), NA),
        med_tdmg = c(median(sm$tdmg, na.rm = TRUE), NA))
}

#' The same statistics from the engine, for one population.
qsp2_engine_summary <- function(traj, pp, co, tri) {
  prep <- qsp2_pop_prep(traj, pp, H = min(traj$H, co$H))
  o <- qsp2_pop_obs(traj, prep, co)
  W <- qsp2_pop_W(traj, tri)
  E <- (W %*% o$S)[1, ]
  c(t2 = E[["t2"]], pos = E[["pos"]], dx = E[["dx"]], dxs = E[["dxs"]], dxc = E[["dxc"]],
    sev = E[["sev"]], sev0 = E[["sev0"]], O1 = E[["o1num"]] / E[["o1den"]],
    O3 = E[["o3num"]] / E[["o3den"]], basep = E[["basep"]], dmg = E[["dmg"]],
    dxpos = E[["dxpos"]], dxpres = E[["dxpres"]], presraised = E[["presraised"]],
    mean_tdx = E[["m1dx"]] / E[["dx"]], med_tdx = qsp2_pop_tq(W, o$M$dx, 0.5)[1],
    med_tfp = qsp2_pop_tq(W, o$M$firstpos, 0.5)[1],
    med_tdmg = qsp2_pop_tq(W, o$M$dmg, 0.5)[1])
}
