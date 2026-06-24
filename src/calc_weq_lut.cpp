/**
 * calc_weq_lut.cpp
 * Phase equilibrium LUT generator for CO2 + lubricant oil (SW68) binary mixture.
 * Computes w_1_eq over T in [283.15, 343.15] K, P in [3.0, 15.0] MPa.
 *
 * Equilibrium equation (vapor = pure CO2, y_1 ≈ 1):
 *   f(x1) = φ1^V(T,P) - x1 * φ1^L(T,P,x1) = 0
 *
 * Solved via bounded secant method. Cross-interaction energy ε_12 uses
 * temperature-dependent k_ij(T) = A·T + B.
 */

#include "saft_mie.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <exception>

// ============================================================================
// Component parameters
// ============================================================================

// -- CO2 (component 0) --
constexpr double MW_CO2   = 44.0;
constexpr double m_CO2    = 1.2625;
constexpr double s_CO2    = 3.5183;    // Angstrom
constexpr double e_CO2    = 261.62;    // K
constexpr double la_CO2   = 5.687;     // lambda_a
constexpr double lr_CO2   = 21.245;    // lambda_r

// -- SW68 lubricant oil (component 1) --
constexpr double MW_OIL   = 634.0;
constexpr double m_OIL    = 12.087;
constexpr double s_OIL    = 4.391;     // Angstrom
constexpr double e_OIL    = 405.77;    // K
constexpr double la_OIL   = 6.0;       // lambda_a
constexpr double lr_OIL   = 18.0;      // lambda_r

// -- Temperature-dependent k_ij coefficients --
// k_ij(T) = A_k * T + B_k
// A_k = -6.719e-7 * MW_oil + 5.339e-4
// B_k =  2.979e-4 * MW_oil - 0.217
constexpr double A_k = -6.719e-7 * MW_OIL + 5.339e-4;
constexpr double B_k =  2.979e-4 * MW_OIL - 0.217;

// ============================================================================
// T-dependent cross-interaction energy
// ============================================================================

double compute_kij(double T) {
    return A_k * T + B_k;
}

/** Recompute e_ij[1] (pair 0,1) and e_ij[2] (pair 1,0) for current T.
 *  Must be called after init_saft_cross_params() every time T changes.
 *  Only e_ij is T-dependent; s_ij, lambda_aij, lambda_rij, c_ij are constant.
 */
void update_eij_for_T(add_args& args, double T) {
    double kij = compute_kij(T);
    double si3  = args.s[0] * args.s[0] * args.s[0];
    double sj3  = args.s[1] * args.s[1] * args.s[1];
    double sij3 = args.s_ij[1] * args.s_ij[1] * args.s_ij[1];
    double eij  = (1.0 - kij) * std::sqrt(si3 * sj3) / sij3
                * std::sqrt(args.e[0] * args.e[1]);
    args.e_ij[1] = eij;   // CO2-oil  pair (0,1)
    args.e_ij[2] = eij;   // oil-CO2  pair (1,0) — symmetric
}

// ============================================================================
// Phase equilibrium residual
// ============================================================================

/** Evaluate f(x1) = φ1^V - x1 * φ1^L(T, P, x1, 1-x1).
 *  Also returns rho_L and phi_1_L for analysis / output.
 */
double equilibrium_residual(double x1, double T, double P,
                            double phi_1_V, add_args& args,
                            double& rho_L_out, double& phi_1_L_out) {
    std::vector<double> x_liq = {x1, 1.0 - x1};
    rho_L_out   = saft_rho_cpp(T, P, x_liq, LIQUID, args);
    phi_1_L_out = saft_fugcoe_cpp(T, rho_L_out, x_liq, args)[0];
    return phi_1_V - x1 * phi_1_L_out;
}

// ============================================================================
// Bounded secant solver for x1
// ============================================================================

struct SecantResult {
    double x1;
    double rho_L;
    double phi_1_L;
    int    iterations;
    bool   converged;
};

SecantResult solve_x1(double T, double P, double phi_1_V,
                      add_args& args, double x_guess) {
    SecantResult res{0.0, 0.0, 0.0, 0, false};
    double rho_L, phi_L;

    // Boundaries for x1
    const double X_LO = 1e-10;
    const double X_HI = 0.9999;

    // Evaluate f at the initial guess
    double x_mid = std::max(X_LO, std::min(X_HI, x_guess));
    double f_mid = equilibrium_residual(x_mid, T, P, phi_1_V, args, rho_L, phi_L);
    if (std::abs(f_mid) < 1e-10) {
        res.x1 = x_mid; res.rho_L = rho_L; res.phi_1_L = phi_L;
        res.iterations = 1; res.converged = true; return res;
    }

    // Determine bracket: find a point with f>0 and a point with f<0
    double x_pos, f_pos, x_neg, f_neg;
    bool have_bracket = false;

    if (f_mid < 0) {
        // f(guess) < 0: use f(X_LO) as positive end (f(X_LO) ≈ φ_V > 0 always)
        double f_lo = equilibrium_residual(X_LO, T, P, phi_1_V, args, rho_L, phi_L);
        if (f_lo > 0) {
            x_pos = X_LO; f_pos = f_lo;
            x_neg = x_mid; f_neg = f_mid;
            have_bracket = true;
        }
    } else {
        // f(guess) > 0: need to find a point where f < 0
        double f_hi = equilibrium_residual(X_HI, T, P, phi_1_V, args, rho_L, phi_L);
        if (std::abs(f_hi) < 1e-10) {
            res.x1 = X_HI; res.rho_L = rho_L; res.phi_1_L = phi_L;
            res.iterations = 2; res.converged = true; return res;
        }
        if (f_hi < 0) {
            x_pos = x_mid; f_pos = f_mid;
            x_neg = X_HI;  f_neg = f_hi;
            have_bracket = true;
        } else {
            // f(X_HI) also > 0 — scan intermediate x1 values
            // Scan BOTH above and below x_mid (f may be non-monotonic
            // when the metastable liquid root exists only at intermediate x1)
            for (double x_test : {0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50}) {
                if (std::abs(x_test - x_mid) < 0.01) continue;  // skip if ~same as x_mid
                double f_test = equilibrium_residual(x_test, T, P, phi_1_V, args, rho_L, phi_L);
                if (std::abs(f_test) < 1e-10) {
                    res.x1 = x_test; res.rho_L = rho_L; res.phi_1_L = phi_L;
                    res.iterations = 3; res.converged = true; return res;
                }
                if (f_test < 0) {
                    // f_mid > 0 (known), f_test < 0 (just checked)
                    // → x_mid is the positive end, x_test is the negative end
                    x_pos = x_mid;  f_pos = f_mid;
                    x_neg = x_test; f_neg = f_test;
                    have_bracket = true;
                    break;
                }
            }
        }
    }

    if (!have_bracket) {
        // No sign change found — full miscibility regime (P >= P_sat)
        res.x1 = X_HI; res.rho_L = rho_L; res.phi_1_L = phi_L;
        double x_ss = phi_1_V / rho_L;  // approximate, not critical
        if (x_ss > X_HI) x_ss = X_HI;
        if (x_ss > x_guess) res.x1 = x_ss;
        res.iterations = 3;
        res.converged = true;
        return res;
    }

    // Secant iteration using bracket [x_pos (f>0), x_neg (f<0)]
    double x_sec0 = x_pos, f_sec0 = f_pos;
    double x_sec1 = x_neg, f_sec1 = f_neg;

    for (int iter = 0; iter < 200; iter++) {
        if (std::abs(f_sec1) < 1e-10 || std::abs(x_sec1 - x_sec0) < 1e-14) {
            res.x1 = x_sec1; res.rho_L = rho_L; res.phi_1_L = phi_L;
            res.iterations = iter + 3; res.converged = true;
            return res;
        }

        double denom = f_sec1 - f_sec0;
        double x2;
        if (std::abs(denom) < 1e-16) {
            x2 = (x_sec0 + x_sec1) * 0.5;
        } else {
            x2 = x_sec1 - f_sec1 * (x_sec1 - x_sec0) / denom;
        }

        // Step-size damping
        double max_step = 0.3 * (x_sec0 + x_sec1);
        if (max_step < 0.01) max_step = 0.01;
        if (x2 < x_sec1 - max_step) x2 = x_sec1 - max_step;
        if (x2 > x_sec1 + max_step) x2 = x_sec1 + max_step;

        // Hard bounds
        if (x2 < X_LO)  x2 = X_LO;
        if (x2 > X_HI)  x2 = X_HI;
        if (!std::isfinite(x2)) x2 = (x_sec0 + x_sec1) * 0.5;

        x_sec0 = x_sec1; f_sec0 = f_sec1;
        x_sec1 = x2;
        f_sec1 = equilibrium_residual(x_sec1, T, P, phi_1_V, args, rho_L, phi_L);
    }

    res.x1 = x_sec1; res.rho_L = rho_L; res.phi_1_L = phi_L;
    res.iterations = 203;
    res.converged = (std::abs(f_sec1) < 1e-8);
    return res;
}

// ============================================================================
// Compute pure-CO2 fugacity coefficient at (T, P)
// ============================================================================

double compute_phi1_V(double T, double P, add_args& args_pure_CO2) {
    double rho_V;
    // Try vapor phase first; fall back to liquid if no vapor root exists
    try {
        rho_V = saft_rho_cpp(T, P, {1.0}, VAPER, args_pure_CO2);
    } catch (...) {
        rho_V = saft_rho_cpp(T, P, {1.0}, LIQUID, args_pure_CO2);
    }
    return saft_fugcoe_cpp(T, rho_V, {1.0}, args_pure_CO2)[0];
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== CO2 + SW68 Phase Equilibrium LUT Generator ===" << std::endl;

    // -- k_ij coefficients --
    std::cout << "k_ij(T) = " << A_k << " * T + " << B_k << std::endl;
    std::cout << "k_ij(283.15 K) = " << compute_kij(283.15) << std::endl;
    std::cout << "k_ij(343.15 K) = " << compute_kij(343.15) << std::endl;

    // ---- Setup binary mixture args ----
    add_args args_mix;
    args_mix.m.push_back(m_CO2);         args_mix.s.push_back(s_CO2);
    args_mix.e.push_back(e_CO2);         args_mix.lambda_a.push_back(la_CO2);
    args_mix.lambda_r.push_back(lr_CO2);
    args_mix.m.push_back(m_OIL);         args_mix.s.push_back(s_OIL);
    args_mix.e.push_back(e_OIL);         args_mix.lambda_a.push_back(la_OIL);
    args_mix.lambda_r.push_back(lr_OIL);
    init_saft_cross_params(args_mix);

    // ---- Setup pure-CO2 args (for vapor-phase fugacity) ----
    add_args args_CO2 = get_single_component(0, args_mix);
    init_saft_cross_params(args_CO2);  // ensure cross-param vectors are sized

    // ---- Grid ----
    constexpr double T_min = 283.15, T_max = 343.15, dT = 1.0;
    constexpr double P_min = 3.0e6,  P_max = 15.0e6,  dP = 0.1e6;

    constexpr int nT = static_cast<int>((T_max - T_min) / dT + 0.5) + 1;
    constexpr int nP = static_cast<int>((P_max - P_min) / dP + 0.5) + 1;
    std::cout << "Grid: " << nT << " x " << nP << " = " << nT * nP << " nodes" << std::endl;

    // ---- Output file ----
    std::ofstream fout("weq_lut.csv");
    fout << std::scientific << std::setprecision(10);
    fout << "T_K,P_Pa,x1,w1_eq,rho_L_molm3,phi1_V,phi1_L,iter,converged" << std::endl;

    // ---- Main double loop ----
    int n_conv = 0, n_notconv = 0, n_skip = 0;
    double x1_prev = 0.3;

    for (int iT = 0; iT < nT; iT++) {
        double T = T_min + iT * dT;

        // Update cross-interaction energy once per isotherm
        update_eij_for_T(args_mix, T);

        // Reset state at start of each isotherm
        x1_prev = 0.3;
        bool entered_full_misc = false;  // once true, stays true for this T

        for (int iP = 0; iP < nP; iP++) {
            double P = P_min + iP * dP;

            // If already in full-miscibility regime, skip solving
            if (entered_full_misc) {
                fout << T << "," << P << ","
                     << "9.9990000000e-01,9.9856102045e-01,"
                     << "nan,nan,nan,0,1" << std::endl;
                continue;
            }

            // --- Pure CO2 fugacity ---
            double phi_1_V;
            try {
                phi_1_V = compute_phi1_V(T, P, args_CO2);
            } catch (const std::exception&) {
                fout << T << "," << P << ",nan,nan,nan,nan,nan,0,0" << std::endl;
                n_skip++;
                continue;
            }

            // --- Solve for x1 ---
            try {
                SecantResult sr = solve_x1(T, P, phi_1_V, args_mix, x1_prev);

                double w1_eq = sr.x1 * MW_CO2 / (sr.x1 * MW_CO2 + (1.0 - sr.x1) * MW_OIL);
                int conv_flag = sr.converged ? 1 : 0;

                fout << T << "," << P << ","
                     << sr.x1 << "," << w1_eq << ","
                     << sr.rho_L << "," << phi_1_V << "," << sr.phi_1_L << ","
                     << sr.iterations << "," << conv_flag << std::endl;

                if (sr.converged) {
                    if (sr.x1 > 0.999) {
                        entered_full_misc = true;  // lock in full miscibility
                    }
                    x1_prev = (sr.x1 > 0.99) ? 0.5 : sr.x1;
                    n_conv++;
                } else {
                    x1_prev = 0.3;  // reset guess for next point
                    n_notconv++;
                }
            } catch (const std::exception&) {
                fout << T << "," << P << ",nan,nan,nan,nan,nan,0,0" << std::endl;
                n_skip++;
            }
        }

        if ((iT + 1) % 10 == 0 || iT == 0 || iT == nT - 1) {
            std::cout << "  T = " << T << " K  (" << iT + 1 << "/" << nT
                      << ")  |  converged: " << n_conv
                      << "  failed: " << n_notconv
                      << "  skipped: " << n_skip << std::endl;
        }
    }

    fout.close();

    std::cout << "\n==============================================" << std::endl;
    std::cout << "  Converged:  " << n_conv    << std::endl;
    std::cout << "  Not converged: " << n_notconv << std::endl;
    std::cout << "  Skipped:    " << n_skip    << std::endl;
    std::cout << "  Output:     weq_lut.csv" << std::endl;
    std::cout << "==============================================" << std::endl;

    return 0;
}
