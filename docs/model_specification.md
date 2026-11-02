# Model specification

This document is Additional file 1 of the article "A quantitative systems pharmacology model of immune checkpoint inhibitor myocarditis does not reproduce clinical onset, incidence and troponin data: a history-matching calibration study" (Abhishek Anil, Aswini Saravanan, Sayantan Shankar Roy and Anand Srinivasan).

It specifies the model as it was calibrated. The code is the definitive implementation:
- `src/qsp2_rhs.c` and `R/qsp2_model.R` hold the equations;
- `R/lib_pop.R` holds the population engine and the observation model;
- `R/lib_cal.R` holds the calibration cohorts and priors.

Parameter names below are those of the code.

## 1. Conventions

- **Time:** days from the first dose. Trajectories are solved to 730 days, and each cohort is observed over its own follow-up.
- **Drug:** amounts in mg and concentrations in µg/mL.
- **Myocardium:** fractions of left-ventricular mass, with mass Mlv = 124 g.
- **Troponin:** ng/L.
- **T cells:** counted as cells.
- **Solver:** LSODA (R package deSolve), relative and absolute tolerances 1e-10 and 1e-12, with the right-hand side compiled from C.

## 2. State variables

| Symbol | Meaning | Unit | Initial value |
|---|---|---|---|
| A1n, A2n | nivolumab, central and peripheral | mg | 0 |
| A1i, A2i | ipilimumab, central and peripheral | mg | 0 |
| A1p, A2p | pembrolizumab, central and peripheral | mg | 0 |
| P | primed autoreactive T cells in lymphoid tissue | cells | 0 |
| E | expanding effector pool | cells | 0 |
| T1 to T_NT | transit chain to the heart (NT compartments) | cells | 0 |
| Eh | effector T cells in the myocardium | cells | 0 |
| Rm | unlicensed resident memory T cells | cells | Rm0 = 1e5 |
| L | licensed resident memory T cells | cells | 0 |
| Mq | activated macrophages | cells | 0 |
| Gr | regulatory T cells in the myocardium | cells | lam_g / del_g |
| X | accumulated exhaustion signal | dimensionless | 0 |
| v | viable myocardium | fraction of LV mass | 1 |
| d | damaged, living myocardium | fraction of LV mass | 0 |
| s | scar | fraction of LV mass | 0 |
| Qs | slow-release troponin pool in damaged tissue | ng | 0 |
| Tn | troponin in the distribution volume | ng/L | Tn0 |

This gives 19 states plus NT transit compartments. Tissue is conserved: v + d + s = 1.

## 3. Module M1: pharmacokinetics

Each antibody j follows a linear two-compartment model with zero-order infusion:

    dA1j/dt = Inf_j(t) − (CL_j/Vc_j) A1j − (Q_j/Vc_j) A1j + (Q_j/Vp_j) A2j
    dA2j/dt = (Q_j/Vc_j) A1j − (Q_j/Vp_j) A2j,        C_j = A1j / Vc_j

Infusions last 0.5 h (nivolumab, pembrolizumab) or 1.5 h (ipilimumab). Clearances scale with (BW/BWref)^aCL and volumes with (BW/BWref)^aV, at a body weight of 80 kg.

| Drug | CL (L/day) | Vc (L) | Q (L/day) | Vp (L) | aCL | aV | BWref (kg) | Source |
|---|---|---|---|---|---|---|---|---|
| nivolumab | 0.2256 | 3.63 | 0.7704 | 2.78 | 0.566 | 0.597 | 80 | Bajaj 2017 |
| ipilimumab | 0.360 | 4.15 | 0.986 | 3.11 | 0.642 | 0.708 | 80 | Feng 2014 |
| pembrolizumab | 0.22 | 3.48 | 0.795 | 4.06 | 0.595 | 0.489 | 76.8 | Ahamadi 2017 |

Nivolumab clearance is held constant, although the source reports a fall over time.

The regimens are:
- nivolumab 240 mg every 2 weeks, which stands for anti-PD-1 monotherapy;
- nivolumab 1 mg/kg plus ipilimumab 3 mg/kg every 3 weeks for four doses, then nivolumab 240 mg every 2 weeks.

## 4. Module M2: receptor occupancy

    Opd1 = (Cn + Cp) / (EC50pd1 + Cn + Cp)
    Octla4 = Ci / (EC50ctla4 + Ci)

## 5. Module M3: priming, expansion and trafficking

S is the patient's autoreactive precursor supply (Section 10).

    prime = lam_p · S · (gC · Octla4 + gPp · Opd1)
    cap   = max(1, S · Emax · (gC · Octla4 + gPp · Opd1))
    dP/dt = prime − (k_act + del_P) P
    dE/dt = k_act P + rho (1 + gP · Opd1) E (1 − E/cap) − (del_E + k_tr) E
    dT1/dt = k_tr (E − T1);   dTk/dt = k_tr (T(k−1) − Tk), k = 2 … NT
    arrive = k_tr · T_NT

- Priming occurs only under checkpoint blockade, so the untreated state is an exact fixed point.
- Ipilimumab alone primes through gC, and an anti-PD-1 antibody alone through gPp.
- Parameter sets must satisfy rho < del_E + k_tr, so that the effector pool contracts when no checkpoint is blocked.

## 6. Module M4: the myocardial immune compartment

    a_f   = (1 + gF · Opd1) / (1 + Gr/Kg)
    lic   = k_lic · Rm · Eh / (Eh + K_lic)
    dEh/dt = phi · arrive + pi_h · Eh · a_f − del_h · Eh − eps_x · X · Eh
    dRm/dt = −lic
    dL/dt  = lic − del_L · L
    dMq/dt = k_m (Eh + L)(1 − Mq/Mmax) − del_m · Mq
    dGr/dt = lam_g + k_g (Eh + L) − del_g · Gr
    eff    = Eh + wL · L + wM · Mq
    kappa  = k_kill · Svul · eff / (Kh + eff) · v
    dX/dt  = k_x · kappa − del_x · X

- Resident memory cells cannot kill until licensed by arriving effector cells, so no injury can precede trafficking.
- Killing saturates with the effector burden, and the exhaustion signal ends the response.
- Myocyte vulnerability Svul is 1.

## 7. Module M5: myocardium and function

    dv/dt = −kappa + r_rep · d
    dd/dt = kappa − (r_rep + r_scar) d
    ds/dt = r_scar · d
    LVEF  = EFmin + (EF0 − EFmin) (v + eta · d)^pexp · (1 − phi_inf · eff/(Kh + eff))

- EFmin = 10%. EF0 is the patient's baseline (Section 10).
- The last factor is a reversible depression of function by inflammation, which saturates where killing saturates.

## 8. Module M6: troponin

    rel    = theta · kappa · Mlv
    dQs/dt = (1 − fc) · rel − k_rel · Qs
    dTn/dt = (fc · rel + k_rel · Qs) / Vd − k_el (Tn − Tn0)

| Quantity | Value | Source or basis |
|---|---|---|
| k_el, troponin T | ln 2 / 134.1 min = 7.44 per day | Kristensen 2024 |
| k_el, troponin I | ln 2 / 239.7 min = 4.16 per day | Kristensen 2024 (assay with the longest half-life) |
| Vd | 5.5 L | Kristensen 2024 (5.1 to 6.2 L) |
| theta | 1.07e5 ng per g destroyed | set so that a Mewton-like injury gives a peak of about 527 ng/L per gram (Mewton 2015) |
| Tn0 | 5 ng/L | nominal reference baseline |

The drug-released troponin T trajectory is computed from the troponin I trajectory by the ratio of clearances. The equations are linear in Tn, so this is exact up to the error checked in the verification report.

## 9. Module M7: tiers, presentation, diagnosis and treatment

- **T1, drug-related injury.** Noise-free troponin, with the reference baseline, exceeds the assay's upper reference limit (URL) at some time.
- **T2, myocarditis.** T1, and myocardial T cells (Eh + L + Rm) exceed f_dx × 7.5e7 cells.
  - 7.5e7 cells is the whole-ventricle equivalent of 7 CD3+ cells/mm² (Caforio 2013).
  - It is converted with an Abercrombie correction: section thickness 4 µm, nuclear diameter 7 µm, left-ventricular volume 118 mL.
  - The fitted factor f_dx absorbs the uncertainty of this conversion and the gap between a histological threshold and a clinical diagnosis.
- **T3, severe myocarditis.** T2 with a nadir LVEF at least 10 percentage points below baseline and below 50%.
- **Clinical presentation.** The first time noise-free troponin exceeds k_sym,i × URL, or LVEF falls below 50%.
- **Diagnosis.** The earlier of two events:
  - the first positive scheduled test at a time when T2 has been reached;
  - clinical presentation, if T2 is reached no later than 7 days after it.

  A positive test before T2 leads to no diagnosis.
- **Treatment.** At diagnosis the checkpoint inhibitor is stopped. Killing halts at t_h = t_dx + L_tx, after which damaged tissue repairs or scars at r_rep and r_scar:

      v(t) = v_h + d_h · r_rep/(r_rep + r_scar) · (1 − exp(−(r_rep + r_scar)(t − t_h)))
      d(t) = d_h · exp(−(r_rep + r_scar)(t − t_h))

  Tests scheduled after a clinical diagnosis are not done. A patient not diagnosed within the horizon is untreated throughout.

## 10. Module M8: patient variation and the virtual population

| Quantity | Distribution |
|---|---|
| Precursor supply S | lognormal, log-scale mean mu, spread sigma |
| Tempo tau | lognormal, median 1, log-scale spread sd_tau |
| Baseline troponin b_i | lognormal, median Tn0, log-scale spread s_b (cohort strata where a source reports them) |
| Symptom threshold k_sym,i | k_sym · exp(s_ksym · z), z standard normal |
| Baseline LVEF EF0 | normal, mean 60%, SD 7% |

- **What tempo does.** Tempo multiplies these rates of the systemic and cardiac response: lam_p, k_act, del_P, rho, del_E, k_tr, pi_h, del_h, eps_x, k_lic, del_L, k_m, del_m, lam_g, k_g, del_g, k_x and del_x. It does not multiply r_rep, r_scar, k_kill or the capacities (K_lic, Kg, Mmax, Kh).
- **Integration over S and tau.** The equations are solved on a grid of S crossed with a grid of tau:
  - S grid: 25 log-spaced nodes from 1e-4 to 1e4 during the search, and 81 nodes for acceptance;
  - tau grid: 25 nodes, 2^k for k = −6 to 6 in steps of 0.5.

  Each cell carries the probability mass between geometric midpoints.
- **Quadrature over the remaining quantities.** Baseline troponin uses 16 slices, placed densely near the cohort's limit. Baseline LVEF and the symptom threshold use 3 points each.
- **Exact population statistics.** Measurement error and background rises enter in closed form, so population statistics are weighted sums with no Monte Carlo error.
- **End-cell rule.**
  - Mass below the lowest S node is allowed only if that node has no T1 at any tempo in either regimen.
  - The highest S cell may hold at most 1%.
  - Tempo end cells hold less than 1% at every allowed sd_tau.
  - A population that breaks the rule, or loses more than 0.1% of solved mass, is not used.
- **Not on the grid.** Myocyte vulnerability, resident-memory seeding, regulatory capacity, body weight and drug clearance are held at their typical values.

## 11. Observation model

- **Measurement.** Tn_obs = (b_i + X(t)) · exp(e), where X(t) = Tn(t) − Tn0 is drug-released troponin and e ~ N(0, 0.09²), independent at each test.
- **Positive result.**
  - A value above the cohort's threshold, or a cohort-specific rule.
  - When the measured baseline is above the limit, a later positive needs at least r_rise × baseline, with r_rise = 1.5 by default.
  - C5 counts every on-treatment value at or above the limit, whatever the baseline.
  - C6 requires a doubling from baseline to at least 28 ng/L.
- **Background rises.**
  - Episodes arrive as a Poisson process at rate lam_bg,i = lam_bg · (b_i / Tn0)^g_bg.
  - Each has a lognormal peak with median m_bg and log-scale spread s_bg = 1 (troponin T scale), and decays at k_el.
  - A scheduled test is positive from background with probability 1 − exp(−lam_bg,i · D(L)), where D(L) = [(ln m_bg − ln L) Φ(z) + s_bg φ(z)] / k_el and z = (ln m_bg − ln L)/s_bg.
- **Assays.**
  - Troponin I assays are related by the ratio of their 99th-percentile URLs to that of the reference assay (Abbott ARCHITECT, 26.2 ng/L). The ratio applies to drug-released troponin, the background threshold and the baseline.
  - This assumes that results at the same multiple of each assay's URL correspond. It is not a measured equivalence of concentrations.
  - Troponin T uses a URL of 14 ng/L.
- **Testing exposure.**
  - The number of tests after the baseline, K, is geometric on 1 to 34 with mean 7.74, from the one source reporting exposure (van den Berg 2024: 2 to 35 measurements, mean 8.74).
  - Test j is done with probability P(K ≥ j) unless a diagnosis came first.
  - All population shares are exact expectations over K.
- **Mixed calendars.** A cohort may be a mixture of calendars, and its statistics are the weighted sums of its parts.

## 12. Calibration cohorts as simulated

| Cohort | Calendar | Assay and threshold | Combination share | Follow-up |
|---|---|---|---|---|
| S0 (C1, C2) | no screening; diagnosis at presentation | hs-cTnI, URL 26.2 ng/L | pure regimens | 365 days |
| S0T (C6 troponin at diagnosis) | no screening | hs-cTnT, URL 14 ng/L | 10% (assumed) | 365 days |
| C3 | monthly, days 30 to 180 | hs-cTnI above 26.2 ng/L | 5% | 180 days |
| C4 | before each infusion: every 21 days (222/298), every 14 or 28 days (38/298 each) | cTnI above 34 ng/L | 7/298 | 730 days, with testing exposure |
| C5 | every 21 days | cTnI at or above 55 ng/L, rule "threshold" | 56/428 | 730 days, with testing exposure |
| C6 (background rises) | every 21 or 28 days (half each) | hs-cTnT, doubling to at least 28 ng/L; baseline median 10.5 ng/L (IQR 6.5 to 15) | 58/164 | 730 days, with testing exposure |

## 13. Parameters

**Fixed from sources.**
- Pharmacokinetics (Section 3).
- Troponin clearance and distribution volume (Section 8).
- Assay URLs and coefficient of variation (Section 11).
- Biopsy threshold and stereological constants (Section 9).
- LVEF cut-off of 50% and fall of 10 points; left-ventricular mass 124 g and volume 118 mL; baseline LVEF distribution.
- Dosing calendars (prescribing information).

**Fixed at nominal values without a direct source.** These are held fixed. The fitted parameters of the same modules (for example lam_p, rho, phi, pi_h and k_kill) can partly absorb them, but they are not identified from data.

| Parameter | Value | Parameter | Value |
|---|---|---|---|
| k_act | 0.2 per day | Mmax | 1e8 cells |
| del_P | 0.05 per day | del_m | 0.1 per day |
| Emax | 1e9 cells | lam_g | 1e3 cells per day |
| del_E | 0.05 per day | k_g | 0.01 per day |
| Kg | 1e5 cells | del_g | 0.05 per day |
| del_L | 0.05 per day | del_x | 0.05 per day |
| k_m | 0.02 per day | Rm0 | 1e5 cells |
| Tn0 | 5 ng/L | s_bg | 1.0 |

**Fitted (36 model parameters).** The priors are listed in Additional file 2 (sheet "Priors") and in `R/lib_cal.R`:

- **Pharmacodynamics and immune response:** lam_p, gC, gPp, rho, gP, gF, k_tr, NT, phi, pi_h, del_h, k_kill, Kh, k_lic, K_lic, wL, wM, eps_x, k_x, EC50pd1, EC50ctla4.
- **Tissue and function:** r_rep, r_scar, eta, pexp, phi_inf.
- **Troponin release:** fc, k_rel.
- **Diagnosis and treatment:** f_dx, k_sym, s_ksym, L_tx.
- **Observation model:** s_b, lam_bg, m_bg, g_bg.

**Population parameters.** Evaluated on a grid for every parameter set:
- mu: −12 to 4 in steps of 0.5;
- sigma: 0.4 to 4.4 in steps of 0.25;
- sd_tau: 0.05, 0.15, 0.3, 0.5, 0.75, 1, 1.25 and 1.5.

## 14. Known simplifications

1. The conversion of the biopsy threshold to cells per left ventricle carries about a factor of three of uncertainty. It is absorbed by f_dx.
2. theta is set by matching a peak troponin per gram, not measured as a release per gram.
3. EC50pd1 has no primary measurement and is fitted with a wide prior.
4. Nivolumab clearance is constant, and body weight is the only pharmacokinetic covariate.
5. Pembrolizumab regimens are represented by nivolumab every 2 weeks. In eight well-fitting parameter sets the two gave identical T2 incidence and timing, but this was not shown for the whole prior.
6. Treatment is idealised as a halt to killing after a lag. Steroid pharmacology and non-response are not modelled.
7. Drug-related troponin trajectories assume dosing continues to the horizon. For a patient whose treatment stops, only testing stops.
8. Testing exposure is geometric, and the assay conversion rests on reference limits. Both are assumptions.
9. The C6 regimen mix is assumed at 10%. In the accepted ensemble, the C6 statistic changes by 1.2% between 0% and 30%.
10. The immune constants in the nominal-value table above have no direct source.
