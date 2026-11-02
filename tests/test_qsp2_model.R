# testthat regression tests for the qsp2 model.
#
# These pin the behaviour that later blocks rely on, so that a change to the model
# cannot pass unnoticed.  Run them with R/run_tests.R; block B4 and every later block
# runs them before doing anything else.

library(testthat)

test_that("the two right-hand sides agree", {
  p <- qsp2_default_pars(); y <- qsp2_init_state(p)
  set.seed(1)
  for (i in 1:50) {
    yr <- y
    yr[] <- abs(yr) + runif(length(yr), 0, 1e4)
    j <- 8L + as.integer(round(p[["NT"]]))
    v <- runif(1, 0, 1); d <- runif(1, 0, 1 - v)
    yr[j + 7] <- v; yr[j + 8] <- d; yr[j + 9] <- 1 - v - d
    # the compiled side returns a named vector, the R side an unnamed one; the test
    # is about the values, so the names are dropped before comparing
    a <- unname(qsp2_rhs_C(0, yr, p)$dy)
    b <- unname(qsp2_rhs_R(0, yr, p)[[1]])
    expect_equal(a, b, tolerance = 1e-12)
  }
})

test_that("the untreated model is a fixed point", {
  p <- qsp2_default_pars(); y <- qsp2_init_state(p)
  r <- qsp2_simulate(p, qsp2_regimen("none"), c(0, 730))
  for (nm in intersect(names(r), names(y)))
    expect_equal(r[[nm]][2], r[[nm]][1], tolerance = 1e-8, info = nm)
})

test_that("priming needs a blocked checkpoint", {
  p <- qsp2_default_pars(); y <- qsp2_init_state(p)
  expect_equal(unname(qsp2_rhs_R(0, y, p)[[1]][7]), 0)
})

test_that("tissue is conserved and scar never falls", {
  p <- qsp2_default_pars(); p[["k_kill"]] <- 0.05; p[["Sus"]] <- 3
  r <- qsp2_simulate(p, qsp2_regimen("combo_nivo1_ipi3"), seq(0, 365, by = 1))
  expect_true(max(abs(r$v + r$d + r$s - 1)) < 1e-9)
  expect_true(all(diff(r$s) >= -1e-12))
})

test_that("no killing means no troponin", {
  p <- qsp2_default_pars(); p[["k_kill"]] <- 0
  r <- qsp2_simulate(p, qsp2_regimen("combo_nivo1_ipi3"), seq(0, 365, by = 5))
  expect_true(max(abs(r$Tn - p[["Tn0"]])) < 1e-9)
})

test_that("ipilimumab alone can still cause injury", {
  p <- qsp2_default_pars(); p[["k_kill"]] <- 0.05; p[["Sus"]] <- 3
  r <- qsp2_simulate(p, qsp2_regimen("ipi_alone"), seq(0, 365, by = 1))
  expect_true(max(r$Tn) > p[["Tn0"]])
  expect_true(max(r$Opd1) < 1e-12)
})

test_that("the published half-lives are reproduced", {
  p <- qsp2_default_pars()
  hl <- function(CL, Vc, Q, Vp) {
    k10 <- CL / Vc; k12 <- Q / Vc; k21 <- Q / Vp
    s <- k10 + k12 + k21; r <- sqrt(s^2 - 4 * k10 * k21)
    log(2) / (0.5 * (s - r))
  }
  expect_equal(hl(p[["CLi"]], p[["Vci"]], p[["Qi"]], p[["Vpi"]]), 14.7, tolerance = 0.05)
  expect_equal(hl(p[["CLp"]], p[["Vcp"]], p[["Qp"]], p[["Vpp"]]), 27.3, tolerance = 0.08)
})

# --- the observation model -------------------------------------------------------------
# Each of these pins a fault that was in the code for the whole first run of this study.

test_that("the assay argument is live, not decorative", {
  p <- qsp2_default_pars(assay = "hs_cTnI")
  p[["Sus"]] <- 20
  reg <- qsp2_regimen("combo_nivo1_ipi3"); tt <- seq(0, 150, by = 1)
  ri <- qsp2_simulate(p, reg, tt, assay = "hs_cTnI")
  rt <- qsp2_simulate(p, reg, tt, assay = "hs_cTnT")
  expect_true(max(abs(ri$Tn - rt$Tn)) > 1e-6 * max(ri$Tn))
  expect_true(max(rt$Tn) < max(ri$Tn))          # troponin T clears about twice as fast
})

test_that("naming no assay leaves the rates in the parameter vector alone", {
  # the B4 dimensional check depends on this: it divides every rate by 24 and would be
  # silently undone if this function put the clearance rate back into days
  p <- qsp2_default_pars(assay = "hs_cTnI")
  p[["k_el"]] <- p[["k_el"]] / 24
  reg <- qsp2_regimen("nivo_q2w")
  r1 <- qsp2_simulate(p, reg, seq(0, 50, by = 1))
  r2 <- qsp2_simulate(p, reg, seq(0, 50, by = 1), assay = "hs_cTnI")
  expect_false(isTRUE(all.equal(max(r1$Tn), max(r2$Tn))))
})

test_that("the myocarditis tier scales with its fitted factor", {
  p <- qsp2_default_pars()
  expect_equal(qsp2_t2_cells(p), QSP2_T2$cells)
  p[["f_dx"]] <- 0.1
  expect_equal(qsp2_t2_cells(p), 0.1 * QSP2_T2$cells)
  q <- qsp2_default_pars(); q <- q[setdiff(names(q), "f_dx")]
  expect_equal(qsp2_t2_cells(q), QSP2_T2$cells)   # absent means unchanged
})

test_that("a cohort carries its own assay and threshold", {
  co <- qsp2_cohort(c(0, 14, 28), "hs_cTnT", 14)
  expect_equal(co$limit, 14)
  expect_equal(co$assay, "hs_cTnT")
  expect_equal(qsp2_cohort(c(0, 14))$limit, QSP2_TROPONIN$hs_cTnI$limit)
  expect_error(qsp2_cohort(c(0, 14), "not_an_assay"))
})

test_that("the chance of a raised result has one definition", {
  # one function holds this arithmetic
  expect_equal(qsp2_p_raised(0.1, 0, 10), 0.1)            # no background, drug only
  expect_equal(qsp2_p_raised(0, 0.01, 10), qsp2_p_background(0.01, 10))
  expect_true(qsp2_p_raised(0.1, 0.01, 25) > qsp2_p_raised(0.1, 0.01, 6))
  expect_true(qsp2_p_raised(0.1, 0.01, 25) <= 1)
  expect_equal(qsp2_p_raised(0.1, 0.01, 4), 1 - 0.9 * 0.99^4)
})

# --- the population engine ---------------------------------------------------------------

test_that("the LVEF criterion needs a fall of 10 points to below 50%", {
  expect_equal(qsp2_g_crit(60), 0.8)
  expect_equal(qsp2_g_crit(70), 40 / 60)                 # below 50% is the binding part
  expect_equal(qsp2_g_crit(52), 1 - 10 / 42)             # the 10-point fall is binding
  expect_equal(qsp2_S_c(0.8), 0)                         # nobody meets it above G = 0.8
  e <- 60 + 7 * qnorm((1:100000 - 0.5) / 100000)
  expect_equal(qsp2_S_c(0.7), mean(qsp2_g_crit(e) > 0.7), tolerance = 1e-4)
})

test_that("severity bands add up to T3, and patients who presented by LVEF are all T3", {
  s <- .sev_c(c(0.3, 0.55, 0.75), 0, 1)
  expect_equal(unname(rowSums(s[, 2:4])), unname(s[, 1]))
  f <- .sev_c(0.9, 0.5, 0.7, forced = TRUE)
  expect_equal(unname(f[1, "t3"]), 1)
})

test_that("the baseline split puts the raised share exactly", {
  co <- qsp2_pop_cohort("t", days = c(21, 42), kin = "TnI", url = 34, limit = 34)
  pp <- list(Tn0 = 5, s_b = 1.1)
  bn <- qsp2_pop_basenodes(co, pp)
  expect_equal(sum(bn$w), 1)
  expect_equal(sum(bn$w[bn$raised]),
               pnorm((log(34) - log(5 * 34 / 26.2)) / 1.1, lower.tail = FALSE))
  expect_true(all(bn$b[bn$raised] > 34) && all(bn$b[!bn$raised] < 34))
})

test_that("background rises fall with the threshold", {
  co <- qsp2_pop_cohort("t", days = c(21, 42), kin = "TnT", url = 14, limit = 14)
  pp <- list(m_bg = 40, s_bg = 1, lam_bg = 0.05)
  p <- qsp2_pop_pbg(c(14, 28, 56), co, pp)
  expect_true(all(diff(p) < 0) && all(p > 0 & p < 1))
})

test_that("population time quantiles put a single day's mass at that day", {
  C <- matrix(cumsum(c(rep(0, 56), 1, rep(0, 10))), nrow = 1)
  expect_equal(qsp2_pop_tq_cum(C, 0.5)[1, 1], 56)
})

test_that("the grid weights of a population sum to one", {
  tr <- list(sgrid = 10^seq(-2, 2, 1), taugrid = c(0.5, 1, 2), ok = rep(TRUE, 15),
             nodes = expand.grid(i_s = 1:5, i_t = 1:3))
  W <- qsp2_pop_W(tr, data.frame(mu = c(0, -1), sigma = c(1, 2), sd_tau = c(0.3, 1)))
  expect_equal(unname(rowSums(W)), c(1, 1))
})

test_that("schedules work for horizons shorter than the quarterly phase and short calendars", {
  # a six-month horizon must not make seq(252, 182, by = 91) fail
  reg <- qsp2_regimen("nivo_q2w", horizon = 182)
  for (s in c("S0", "S1", "S2", "S3", "S4")) expect_true(all(qsp2_schedule(s, reg, 182) <= 182))
  expect_equal(max(qsp2_schedule("S4", reg, 182)), 168)
  expect_equal(qsp2_schedule("S1", qsp2_regimen("ipi_alone"), 365), c(0, 21, 42, 63))
})

test_that("inflammation depresses function reversibly", {
  p <- qsp2_default_pars()
  expect_equal(qsp2_phi_inf(p), 0.4)
  q <- p[setdiff(names(p), "phi_inf")]
  expect_equal(qsp2_phi_inf(q), 0)                       # absent means zero
  # with no burden the map is the tissue-only one; with burden it is lower, and it recovers
  expect_equal(qsp2_lvef(1, 0, p, sat = 0), QSP2_MYO$EF0)
  expect_lt(qsp2_lvef(1, 0, p, sat = 1), qsp2_lvef(1, 0, p, sat = 0))
  expect_equal(qsp2_lvef(1, 0, p, sat = 1), 10 + 50 * (1 - 0.4))
  expect_equal(qsp2_lvef(0.9, 0.05, q, sat = 1), qsp2_lvef(0.9, 0.05, q, sat = 0))
})

test_that("the effector burden that depresses function is the one that drives killing", {
  p <- qsp2_default_pars(); p[["Sus"]] <- 5
  r <- qsp2_simulate(p, qsp2_regimen("combo_nivo1_ipi3"), seq(0, 120, by = 4), assay = "hs_cTnI")
  s <- qsp2_inflam_sat(r, p)
  expect_true(all(s >= 0 & s <= 1))
  eff <- r$Eh + p[["wL"]] * r$L + p[["wM"]] * r$Mq
  expect_equal(s, eff / (p[["Kh"]] + eff))
  expect_equal(s[1], 0)                                   # nothing in the heart before treatment
})

test_that("dosing calendars work for horizons shorter than a later phase", {
  # found by a block B4 test of a finite induction course: seq(84, 63, by = 14) counted backwards
  for (k in c("nivo_q2w", "nivo_q4w", "pembro_q3w", "pembro_q6w", "combo_nivo1_ipi3",
              "combo_nivo3_ipi1", "ipi_alone")) {
    for (h in c(7, 63, 84, 365)) {
      r <- qsp2_regimen(k, horizon = h)
      expect_true(all(r$day <= h), info = paste(k, h))
    }
  }
  expect_equal(nrow(qsp2_regimen("combo_nivo1_ipi3", horizon = 63)), 8L)   # 4 doses of each drug
  expect_equal(max(qsp2_regimen("ipi_alone", horizon = 30)$day), 21)
})

test_that("the peak troponin histogram is at or above the one at presentation", {
  # C6 reports each patient's highest troponin; the engine keeps both the highest value and
  # the value at diagnosis, and they must not be the same object
  p <- qsp2_default_pars(assay = "hs_cTnI"); p[["Sus"]] <- 5; p[["k_kill"]] <- 1e-3
  tr <- qsp2_pop_traj(p, qsp2_regimen("combo_nivo1_ipi3"), sgrid = c(1, 10), taugrid = 1)
  pp <- qsp2_pop_pp(p)
  o <- qsp2_pop_obs(tr, qsp2_pop_prep(tr, pp), qsp2_pop_cohort("S0T", kin = "TnT", url = 14))
  expect_true(all(dim(o$trop) == dim(o$troppk)))
  med <- function(h) { w <- which(cumsum(h) >= 0.5 * sum(h))[1]; QSP2_TROP_BINS[w] }
  for (k in which(rowSums(o$trop) > 0))
    expect_gte(med(o$troppk[k, ]), med(o$trop[k, ]))
})

test_that("background rises follow the patient's own baseline", {
  co <- qsp2_pop_cohort("t", days = c(21, 42), kin = "TnT", url = 14, limit = 14)
  pp <- list(m_bg = 40, s_bg = 1, lam_bg = 0.05, Tn0 = 5, g_bg = 1)
  ref <- pp$Tn0 * co$url / QSP2_REF_URL             # the reference baseline on this assay
  # a sicker patient has more background; at g_bg = 0 everyone is the same
  expect_gt(qsp2_pop_pbg(14, co, pp, b = 4 * ref), qsp2_pop_pbg(14, co, pp, b = ref))
  expect_lt(qsp2_pop_pbg(14, co, pp, b = 0.25 * ref), qsp2_pop_pbg(14, co, pp, b = ref))
  expect_equal(qsp2_pop_pbg(14, co, pp, b = ref), qsp2_pop_pbg(14, co, pp))
  pp0 <- pp; pp0$g_bg <- 0
  expect_equal(qsp2_pop_pbg(14, co, pp0, b = 9 * ref), qsp2_pop_pbg(14, co, pp0))
  pp$g_bg <- NULL                                    # absent behaves as zero
  expect_equal(qsp2_pop_pbg(14, co, pp, b = 9 * ref), qsp2_pop_pbg(14, co, pp))
})

test_that("the persisting-depression alternative lies between treatment and no treatment", {
  # a course in which the depression does not resolve with the infiltrate can never be better
  # than the resolving one nor worse than none
  p <- qsp2_default_pars(assay = "hs_cTnI")
  for (nm in names(m <- list(f_dx = 0.1, k_kill = 1e-3, eta = 0.3, pexp = 2.5, phi_inf = 0.5)))
    p[[nm]] <- m[[nm]]
  tr <- qsp2_pop_traj(p, qsp2_regimen("combo_nivo1_ipi3"), sgrid = c(3, 30), taugrid = 1)
  pp <- qsp2_pop_pp(p)
  a <- qsp2_pop_prep(tr, pp, H = 365, persist = FALSE)$gpost
  b <- qsp2_pop_prep(tr, pp, H = 365, persist = TRUE)$gpost
  k <- which(tr$ok)
  expect_true(all(b[k, ] <= a[k, ] + 1e-12))
  expect_true(all(b[k, ] >= tr$cumG[k, 366] - 1e-12))
  expect_true(any(abs(a[k, ] - b[k, ]) > 1e-6))
})

test_that("tempo scales the clone's response but not the tissue", {
  p <- qsp2_default_pars()
  q <- qsp2_apply_tempo(p, 3)
  for (nm in QSP2_TAU_RATES) expect_equal(q[[nm]], 3 * p[[nm]], info = nm)
  # the myocardium, the killing constant and the capacities keep their own values
  for (nm in c("r_rep", "r_scar", "k_kill", "K_lic", "Kg", "Mmax", "Kh", "theta", "k_rel", "k_el"))
    expect_equal(q[[nm]], p[[nm]], info = nm)
  # the resting state is a ratio of rates that carry the same factor, so it does not move
  expect_equal(qsp2_init_state(q), qsp2_init_state(p), tolerance = 1e-12)
  # and the cardiac compartment is included
  expect_true(all(c("del_h", "pi_h", "k_lic", "del_g") %in% QSP2_TAU_RATES))
})

# --- assays, exposure, end cells, C1 and the search -----------------------------------------

.t18 <- local({
  p <- qsp2_default_pars(assay = "hs_cTnI"); p[["Sus"]] <- 5; p[["k_kill"]] <- 1e-3
  p[["lam_bg"]] <- 0.02
  tr <- qsp2_pop_traj(p, qsp2_regimen("combo_nivo1_ipi3"), sgrid = c(0.3, 3, 30), taugrid = c(0.5, 1))
  pp <- qsp2_pop_pp(p)
  list(tr = tr, pp = pp, prep = qsp2_pop_prep(tr, pp))
})
.obs18 <- function(...) qsp2_pop_obs(.t18$tr, .t18$prep, qsp2_pop_cohort("t", days = seq(21, 365, 21), ...))

test_that("troponin I assays are related by the ratio of their limits", {
  a <- .obs18(kin = "TnI", url = 26.2, limit = 26.2)
  b <- .obs18(kin = "TnI", url = 55, limit = 55)
  # at the same multiple of each assay's own limit, nothing changes
  expect_equal(a$S[, "pos"], b$S[, "pos"], tolerance = 1e-12)
  expect_equal(a$S[, "dx"], b$S[, "dx"], tolerance = 1e-12)
  # on one assay, a higher decision threshold finds fewer positives: the limit is read
  h <- .obs18(kin = "TnI", url = 26.2, limit = 52.4)
  expect_true(all(h$S[, "pos"] <= a$S[, "pos"] + 1e-12) && any(h$S[, "pos"] < a$S[, "pos"] - 1e-6))
  # with no conversion, 55 ng/L is twice the reference limit and finds fewer
  n <- .obs18(kin = "TnI", url = 55, limit = 55, conv = "none")
  expect_true(any(n$S[, "pos"] < a$S[, "pos"] - 1e-6))
})

test_that("the testing exposure has the reported mean and keeps only tests within follow-up", {
  pk <- qsp2_pop_expo_pk(list(mean = 7.74, max = 34), 40)
  expect_equal(sum(pk), 1, tolerance = 1e-12)
  expect_equal(pk[1], 0)                                   # every patient has a test after baseline
  expect_equal(sum((0:40) * pk), 7.74, tolerance = 1e-8)
  short <- qsp2_pop_expo_pk(list(mean = 7.74, max = 34), 5)
  expect_equal(sum(short), 1, tolerance = 1e-12)           # the rest is lumped on the last test
  expect_equal(qsp2_pop_expo_pk(NULL, 6), c(rep(0, 6), 1))
})

test_that("a shorter exposure means fewer tests and never more positives", {
  f <- .obs18(kin = "TnI", url = 34, limit = 34)
  e <- .obs18(kin = "TnI", url = 34, limit = 34, expo = list(mean = 3, max = 34))
  # never more tests, and fewer wherever patients are not all diagnosed by the first test
  expect_true(all(e$S[, "ntests"] <= f$S[, "ntests"] + 1e-12) &&
                any(e$S[, "ntests"] < f$S[, "ntests"] - 1e-3))
  expect_true(all(e$S[, "pos"] <= f$S[, "pos"] + 1e-12))
  expect_equal(rowSums(e$M$firstpos), e$S[, "pos"], tolerance = 1e-10)
  expect_equal(rowSums(e$M$dx), e$S[, "dx"], tolerance = 1e-10)
})

test_that("a mixture of calendars is the weighted sum of its parts", {
  c1 <- qsp2_pop_cohort("a", days = seq(14, 365, 14), kin = "TnI", url = 34, limit = 34)
  c2 <- qsp2_pop_cohort("b", days = seq(28, 365, 28), kin = "TnI", url = 34, limit = 34)
  m <- qsp2_pop_obs_any(.t18$tr, .t18$prep, qsp2_pop_cohort_mix("m", list(c1, c2), c(0.3, 0.7)))
  o1 <- qsp2_pop_obs(.t18$tr, .t18$prep, c1); o2 <- qsp2_pop_obs(.t18$tr, .t18$prep, c2)
  expect_equal(m$S, 0.3 * o1$S + 0.7 * o2$S, tolerance = 1e-12)
  expect_equal(m$M$firstpos, 0.3 * o1$M$firstpos + 0.7 * o2$M$firstpos, tolerance = 1e-12)
})

test_that("R1 is the latent phase alone", {
  st <- data.frame(R1_p7_mono = c(0.96, 0.96, 0.90), R1_med_mono = c(12, 60, 60))
  expect_equal(qsp2_req_R1(st), c(TRUE, TRUE, FALSE))
})

test_that("C1 uses cases with both dates complete to the day", {
  c1 <- qsp2_c1_bootstrap(B = 10)
  expect_equal(c(c1$n_combo, c1$n_mono), c(262, 988))
  expect_equal(c1$q_combo, c(14.25, 20, 39.75)); expect_equal(c1$q_mono, c(20, 30, 89))
  old <- qsp2_c1_bootstrap(B = 10, complete = FALSE)          # completed dates kept as well
  expect_equal(c(old$n_combo, old$n_mono), c(281, 1108))
  expect_equal(old$q_combo, c(14, 21, 41)); expect_equal(old$q_mono, c(19, 31, 94.5))
  tg <- qsp2_cal_targets(c1)
  expect_equal(tg$obs[tg$id == "C1_q75_mono"], 89)            # read from the data, not typed
  expect_equal(nrow(tg), 26L)
  expect_false("C5_p_presented" %in% tg$id)                   # a sensitivity analysis
  expect_equal(qsp2_cal_targets_sens()$id, "C5_p_presented")
})

test_that("each point gets its own box in the search", {
  set.seed(3)
  pr <- list(a = c(0, 1, "unif"), b = c(1, 1000, "log"))
  keep <- data.table(a = c(0.1, 0.9), b = c(2, 500))
  d <- qsp2_lhs_near(400, pr, keep, 0.05)
  near1 <- abs(d$a - 0.1) <= 0.05 + 1e-12 & abs(log(d$b / 2)) <= 0.05 * log(1000) + 1e-12
  near2 <- abs(d$a - 0.9) <= 0.05 + 1e-12 & abs(log(d$b / 500)) <= 0.05 * log(1000) + 1e-12
  expect_true(all(near1 | near2))                             # nothing in between the two
  expect_equal(sum(near1), 200L); expect_equal(sum(near2), 200L)
})

test_that("population mass beyond the S grid is allowed only where it is harmless", {
  mk <- function(t1) list(t1ref = t1, nodes = data.frame(i_s = 1:3), sgrid = c(1e-4, 1, 1e4))
  tri <- data.frame(mu = c(-12, 0, 6), sigma = c(3, 1, 1), sd_tau = 0.5)
  free <- qsp2_cal_edges(list(tr = list(mono = mk(c(NA, 5, 3)), combo = mk(c(NA, 4, 2)))), tri)
  expect_equal(unname(free[, "edge_ok"]), c(1, 1, 0))         # low tail fine, high tail not
  hurt <- qsp2_cal_edges(list(tr = list(mono = mk(c(NA, 5, 3)), combo = mk(c(40, 4, 2)))), tri)
  expect_equal(unname(hurt[, "edge_ok"]), c(0, 1, 0))         # the lowest node has injury
})
