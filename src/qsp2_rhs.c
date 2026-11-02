/* ---------------------------------------------------------------------------
 * qsp2_rhs.c  -  right-hand side of the qsp2 model.
 *
 * This file implements, term for term, docs/model_specification.md.  Section
 * numbers in the comments refer to that document.
 * An independent implementation of the same equations is written in R in
 * R/qsp2_rhs_reference.R; block B4 requires the two to agree to 1e-10, which is
 * the check that this code matches the specification.
 *
 * States (Section 2), in the order used throughout:
 *    0  A1n   nivolumab central          mg
 *    1  A2n   nivolumab peripheral       mg
 *    2  A1i   ipilimumab central         mg
 *    3  A2i   ipilimumab peripheral      mg
 *    4  A1p   pembrolizumab central      mg
 *    5  A2p   pembrolizumab peripheral   mg
 *    6  P     primed cells               cells
 *    7  E     expanding effector pool    cells
 *    8..8+NT-1   transit chain T1..TNT   cells
 *    +0   Eh    cardiac effectors        cells
 *    +1   Rm    resident memory, unlicensed
 *    +2   L     licensed resident memory
 *    +3   Mq    activated macrophages
 *    +4   Gr    regulatory T cells
 *    +5   X     exhaustion signal        dimensionless
 *    +6   v     viable myocardium        fraction
 *    +7   d     damaged myocardium       fraction
 *    +8   s     scar                     fraction
 *    +9   Qs    slow troponin pool       ng
 *    +10  Tn    troponin in blood        ng/L
 *
 * The transit chain length NT is a parameter, so the state count is 19 + NT.
 * Infusions are handled outside the right-hand side, by integrating piecewise
 * between dose boundaries and adding each dose as a zero-order input; the rate
 * for the interval is passed in as inf[] (mg/day into the central compartment).
 * ---------------------------------------------------------------------------
 */

#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* ---- parameters, filled by initmod from the deSolve parameter vector ---- */
#define NPAR 44
static double p[NPAR];
static int    NT;          /* transit chain length, p[43] rounded */

/* index the parameter vector by name, in the order R passes them */
#define CLn   p[0]
#define Vcn   p[1]
#define Qn    p[2]
#define Vpn   p[3]
#define CLi   p[4]
#define Vci   p[5]
#define Qi    p[6]
#define Vpi   p[7]
#define CLp   p[8]
#define Vcp   p[9]
#define Qp    p[10]
#define Vpp   p[11]
#define EC50pd1   p[12]
#define EC50ctla4 p[13]
#define lam_p p[14]
#define gC    p[15]
#define k_act p[16]
#define del_P p[17]
#define rho   p[18]
#define gP    p[19]
#define Emax  p[20]
#define del_E p[21]
#define k_tr  p[22]
#define phi   p[23]
#define pi_h  p[24]
#define gF    p[25]
#define Kg    p[26]
#define del_h p[27]
#define eps_x p[28]
#define k_lic p[29]
#define K_lic p[30]
#define del_L p[31]
#define k_m   p[32]
#define Mmax  p[33]
#define del_m p[34]
#define lam_g p[35]
#define k_g   p[36]
#define del_g p[37]
#define k_x   p[38]
#define del_x p[39]
#define k_kill p[40]
#define Kh    p[41]
#define wL    p[42]
/* p[43] is NT; the remaining model parameters arrive in a second block */

#define NPAR2 16
static double q[NPAR2];
#define wM     q[0]
#define Sus    q[1]
#define r_rep  q[2]
#define r_scar q[3]
#define theta  q[4]
#define fc     q[5]
#define k_rel  q[6]
#define k_el   q[7]
#define Vd     q[8]
#define Mlv    q[9]
#define Tn0    q[10]
#define infn   q[11]   /* nivolumab infusion rate, mg/day */
#define infi   q[12]   /* ipilimumab */
#define infp   q[13]   /* pembrolizumab */
#define gPp    q[14]   /* effect of PD-1 blockade on priming */
#define Svul   q[15]   /* myocyte vulnerability */

void qsp2_init(void (*odeparms)(int *, double *))
{
    int n = NPAR + NPAR2;
    double all[NPAR + NPAR2];
    odeparms(&n, all);
    for (int k = 0; k < NPAR;  k++) p[k] = all[k];
    for (int k = 0; k < NPAR2; k++) q[k] = all[NPAR + k];
    NT = (int) (p[43] + 0.5);
    if (NT < 1) NT = 1;
}

/* self-test kept from block B0: proves the .C interface and the build work */
void qsp2_selftest(double *x, int *n, double *y)
{
    for (int i = 0; i < *n; i++) y[i] = 2.0 * x[i] + 1.0;
}

void qsp2_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    const int iT  = 8;          /* first transit compartment */
    const int iEh = 8 + NT;
    const int iRm = iEh + 1, iL = iEh + 2, iMq = iEh + 3, iGr = iEh + 4, iX = iEh + 5;
    const int iv  = iEh + 6, id = iEh + 7, is = iEh + 8, iQs = iEh + 9, iTn = iEh + 10;

    /* ---- M1 pharmacokinetics (Section 3) ---------------------------------- */
    double Cn = y[0] / Vcn, Ci = y[2] / Vci, Cp = y[4] / Vcp;

    ydot[0] = infn - (CLn / Vcn) * y[0] - (Qn / Vcn) * y[0] + (Qn / Vpn) * y[1];
    ydot[1] =               (Qn / Vcn) * y[0] - (Qn / Vpn) * y[1];
    ydot[2] = infi - (CLi / Vci) * y[2] - (Qi / Vci) * y[2] + (Qi / Vpi) * y[3];
    ydot[3] =               (Qi / Vci) * y[2] - (Qi / Vpi) * y[3];
    ydot[4] = infp - (CLp / Vcp) * y[4] - (Qp / Vcp) * y[4] + (Qp / Vpp) * y[5];
    ydot[5] =               (Qp / Vcp) * y[4] - (Qp / Vpp) * y[5];

    /* ---- M2 receptor occupancy (Section 4) -------------------------------- */
    double Cpd1  = Cn + Cp;
    double Opd1  = Cpd1 / (EC50pd1 + Cpd1);
    double Octla = Ci   / (EC50ctla4 + Ci);

    /* ---- M3 priming, expansion, transit (Section 5) ------------------------ */
    double P = y[6], E = y[7];
    /* Priming happens only under blockade, through either receptor, so the
     * untreated model is an exact fixed point.                                        */
    /* Sus is the patient's autoreactive precursor supply, so it drives
     * the whole cascade and therefore the size of the infiltrate as well as the
     * damage.  Myocyte vulnerability is a separate axis, Svul.                   */
    double prime = lam_p * Sus * (gC * Octla + gPp * Opd1);
    double grow  = rho * (1.0 + gP * Opd1);

    ydot[6] = prime - (k_act + del_P) * P;
    /* The expansion ceiling is the patient's own and depends on the drug.  A ceiling
     * shared by everyone would not let the infiltrate vary between patients, and a
     * ceiling independent of the drug would bring both regimens to the same final
     * infiltrate, so the incidence ratio between them would be exactly one whatever
     * the parameters.  The ceiling is proportional to the same priming drive, so
     * blocking either receptor raises it and the data decide by how much.          */
    double cap = Sus * Emax * (gC * Octla + gPp * Opd1);
    if (cap < 1.0) cap = 1.0;          /* with no drug the pool is empty anyway */
    ydot[7] = k_act * P + grow * E * (1.0 - E / cap) - (del_E + k_tr) * E;

    ydot[iT] = k_tr * (E - y[iT]);
    for (int k = 1; k < NT; k++) ydot[iT + k] = k_tr * (y[iT + k - 1] - y[iT + k]);
    double arrive = k_tr * y[iT + NT - 1];

    /* ---- M4 the myocardium's immune compartment (Section 6) ---------------- */
    double Eh = y[iEh], Rm = y[iRm], L = y[iL], Mq = y[iMq], Gr = y[iGr], X = y[iX];
    double a_f = (1.0 + gF * Opd1) / (1.0 + Gr / Kg);

    double lic = k_lic * Rm * (Eh / (Eh + K_lic));

    ydot[iEh] = phi * arrive + pi_h * Eh * a_f - del_h * Eh - eps_x * X * Eh;
    ydot[iRm] = -lic;
    ydot[iL]  =  lic - del_L * L;
    ydot[iMq] = k_m * (Eh + L) * (1.0 - Mq / Mmax) - del_m * Mq;
    ydot[iGr] = lam_g + k_g * (Eh + L) - del_g * Gr;

    /* ---- killing (Section 6) and M5 myocardium (Section 7) ------------------ */
    double v = y[iv], d = y[id];
    double eff   = Eh + wL * L + wM * Mq;
    double kappa = k_kill * Svul * (eff / (Kh + eff)) * v;

    ydot[iX] = k_x * kappa - del_x * X;

    ydot[iv] = -kappa + r_rep * d;
    ydot[id] =  kappa - (r_rep + r_scar) * d;
    ydot[is] =  r_scar * d;

    /* ---- M6 troponin (Section 8) ------------------------------------------- */
    /* kappa is a fraction of left-ventricular mass per day, so kappa * Mlv is
     * grams per day and theta * kappa * Mlv is ng per day released.           */
    double rel = theta * kappa * Mlv;
    ydot[iQs] = (1.0 - fc) * rel - k_rel * y[iQs];
    ydot[iTn] = (fc * rel + k_rel * y[iQs]) / Vd - k_el * (y[iTn] - Tn0);

    /* ---- reported quantities ------------------------------------------------ */
    yout[0] = Cn;
    yout[1] = Ci;
    yout[2] = Cp;
    yout[3] = Opd1;
    yout[4] = Octla;
    yout[5] = kappa;
    yout[6] = Eh + L + Rm;      /* total cardiac T cells, for the tier T2 test */
    yout[7] = v + d + y[is];    /* tissue total, so conservation is visible in the output */
}
