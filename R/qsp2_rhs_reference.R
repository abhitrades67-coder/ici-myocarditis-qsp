# ---------------------------------------------------------------------------
# qsp2_rhs_reference.R  -  an independent right-hand side, written in R.
#
# This is deliberately a SECOND implementation of docs/model_specification.md,
# written from the specification rather than from src/qsp2_rhs.c.  Block B4
# requires the two to agree to 1e-10 at a thousand random states.  If they
# disagree, at least one of them does not match the specification, which is the
# error this check exists to find.  It is not used in any analysis; it exists to
# test the compiled code.
#
# The equation numbering follows the specification, not the C file, and the terms
# are written out in the order the specification writes them.
# ---------------------------------------------------------------------------

qsp2_rhs_R <- function(t, y, p) {
  NT <- as.integer(round(p[["NT"]]))
  iT  <- 9L                      # R is 1-based: states 1..8 are PK, P and E
  iEh <- 8L + NT + 1L
  iRm <- iEh + 1L; iL <- iEh + 2L; iMq <- iEh + 3L; iGr <- iEh + 4L; iX <- iEh + 5L
  iv  <- iEh + 6L; id <- iEh + 7L; is <- iEh + 8L; iQs <- iEh + 9L; iTn <- iEh + 10L

  dy <- numeric(length(y))

  # Section 3: pharmacokinetics, three independent two-compartment models
  Cn <- y[1] / p[["Vcn"]]; Ci <- y[3] / p[["Vci"]]; Cp <- y[5] / p[["Vcp"]]
  dy[1] <- p[["infn"]] - (p[["CLn"]] / p[["Vcn"]]) * y[1] -
           (p[["Qn"]] / p[["Vcn"]]) * y[1] + (p[["Qn"]] / p[["Vpn"]]) * y[2]
  dy[2] <- (p[["Qn"]] / p[["Vcn"]]) * y[1] - (p[["Qn"]] / p[["Vpn"]]) * y[2]
  dy[3] <- p[["infi"]] - (p[["CLi"]] / p[["Vci"]]) * y[3] -
           (p[["Qi"]] / p[["Vci"]]) * y[3] + (p[["Qi"]] / p[["Vpi"]]) * y[4]
  dy[4] <- (p[["Qi"]] / p[["Vci"]]) * y[3] - (p[["Qi"]] / p[["Vpi"]]) * y[4]
  dy[5] <- p[["infp"]] - (p[["CLp"]] / p[["Vcp"]]) * y[5] -
           (p[["Qp"]] / p[["Vcp"]]) * y[5] + (p[["Qp"]] / p[["Vpp"]]) * y[6]
  dy[6] <- (p[["Qp"]] / p[["Vcp"]]) * y[5] - (p[["Qp"]] / p[["Vpp"]]) * y[6]

  # Section 4: occupancy.  Nivolumab and pembrolizumab compete for PD-1.
  Cpd1  <- Cn + Cp
  Opd1  <- Cpd1 / (p[["EC50pd1"]] + Cpd1)
  Octla <- Ci / (p[["EC50ctla4"]] + Ci)

  # Section 5: priming, expansion, transit
  P <- y[7]; E <- y[8]
  # no constitutive priming; either receptor can permit it
  # susceptibility is the precursor supply, so it scales priming
  dy[7] <- p[["lam_p"]] * p[["Sus"]] * (p[["gC"]] * Octla + p[["gPp"]] * Opd1) -
           (p[["k_act"]] + p[["del_P"]]) * P
  # the ceiling is the patient's own and rises with blockade
  cap <- max(1, p[["Sus"]] * p[["Emax"]] * (p[["gC"]] * Octla + p[["gPp"]] * Opd1))
  dy[8] <- p[["k_act"]] * P +
           p[["rho"]] * (1 + p[["gP"]] * Opd1) * E * (1 - E / cap) -
           (p[["del_E"]] + p[["k_tr"]]) * E
  dy[iT] <- p[["k_tr"]] * (E - y[iT])
  if (NT > 1L) for (k in seq_len(NT - 1L))
    dy[iT + k] <- p[["k_tr"]] * (y[iT + k - 1L] - y[iT + k])
  arrive <- p[["k_tr"]] * y[iT + NT - 1L]

  # Section 6: the myocardium's immune compartment
  Eh <- y[iEh]; Rm <- y[iRm]; L <- y[iL]; Mq <- y[iMq]; Gr <- y[iGr]; X <- y[iX]
  a_f <- (1 + p[["gF"]] * Opd1) / (1 + Gr / p[["Kg"]])
  lic <- p[["k_lic"]] * Rm * (Eh / (Eh + p[["K_lic"]]))

  dy[iEh] <- p[["phi"]] * arrive + p[["pi_h"]] * Eh * a_f -
             p[["del_h"]] * Eh - p[["eps_x"]] * X * Eh
  dy[iRm] <- -lic
  dy[iL]  <-  lic - p[["del_L"]] * L
  dy[iMq] <- p[["k_m"]] * (Eh + L) * (1 - Mq / p[["Mmax"]]) - p[["del_m"]] * Mq
  dy[iGr] <- p[["lam_g"]] + p[["k_g"]] * (Eh + L) - p[["del_g"]] * Gr

  # Section 6 killing and Section 7 myocardium
  v <- y[iv]; d <- y[id]
  eff <- Eh + p[["wL"]] * L + p[["wM"]] * Mq
  kappa <- p[["k_kill"]] * p[["Svul"]] * (eff / (p[["Kh"]] + eff)) * v

  dy[iX] <- p[["k_x"]] * kappa - p[["del_x"]] * X
  dy[iv] <- -kappa + p[["r_rep"]] * d
  dy[id] <-  kappa - (p[["r_rep"]] + p[["r_scar"]]) * d
  dy[is] <-  p[["r_scar"]] * d

  # Section 8: troponin, fast blood compartment plus a slow release pool
  rel <- p[["theta"]] * kappa * p[["Mlv"]]
  dy[iQs] <- (1 - p[["fc"]]) * rel - p[["k_rel"]] * y[iQs]
  dy[iTn] <- (p[["fc"]] * rel + p[["k_rel"]] * y[iQs]) / p[["Vd"]] -
             p[["k_el"]] * (y[iTn] - p[["Tn0"]])

  list(dy, c(Cn = Cn, Ci = Ci, Cp = Cp, Opd1 = Opd1, Octla4 = Octla,
             kappa = kappa, Tcells_heart = Eh + L + Rm, tissue_total = v + d + y[is]))
}

#' Call the compiled right-hand side once, at a single state, for comparison with
#' qsp2_rhs_R.  deSolve::DLLfunc evaluates the derivatives without integrating, so
#' this is a direct term-by-term comparison and not a comparison of two solvers.
qsp2_rhs_C <- function(t, y, p) {
  load_qsp2_dll()
  deSolve::DLLfunc(func = "qsp2_derivs", times = t, y = y, parms = p[QSP2_PAR_NAMES],
                   dllname = "qsp2_rhs", initfunc = "qsp2_init",
                   nout = QSP2_NOUT, outnames = QSP2_OUTNAMES)
}
