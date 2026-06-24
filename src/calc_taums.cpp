/**
 * calc_taums.cpp
 * Step 2: Mesoscopic transport constitutive calculation.
 * Computes τ_ms (characteristic mass transfer time) over the (T,P,w1) grid.
 *
 * Theory chain:
 *   1. SAFT EOS        → ρ_m (mixture density)
 *   2. Free Volume Th.  → D₁₂ (dynamic diffusion coefficient)
 *   3. Viscosity model  → μ_m (mixture viscosity)
 *   4. Dimensionless #  → Re, Sc, Sh
 *   5. τ_ms = L² / (D₁₂ · Sh)
 *
 * All physical constants are at the top of the file with literature sources
 * noted. Replace with experimentally-fitted values when available.
 */

#include "saft_mie.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <exception>
#include <algorithm>

// ============================================================================
// Component SAFT parameters  (same as Step 1)
// ============================================================================
constexpr double MW_CO2  = 44.0;
constexpr double m_CO2   = 1.2625;
constexpr double s_CO2   = 3.5183;
constexpr double e_CO2   = 261.62;
constexpr double la_CO2  = 5.687;
constexpr double lr_CO2  = 21.245;

constexpr double MW_OIL  = 634.0;
constexpr double m_OIL   = 12.087;
constexpr double s_OIL   = 4.391;
constexpr double e_OIL   = 405.77;
constexpr double la_OIL  = 6.0;
constexpr double lr_OIL  = 18.0;

// T-dependent k_ij
constexpr double A_k = -6.719e-7 * MW_OIL + 5.339e-4;
constexpr double B_k =  2.979e-4 * MW_OIL - 0.217;

// ============================================================================
// Free Volume Theory parameters  [PLACEHOLDERS — replace with experiment]
// ============================================================================

// D₀: pre-exponential factor for CO₂ diffusion in heavy oil
// Typical range: 1e-9 to 1e-7 m²/s
// Ref: Vrentas & Duda, J. Polym. Sci. (1977); Areerat et al. (2002)
constexpr double D0_CO2 = 5.0e-8;    // m²/s

// γ: overlap factor accounting for same free volume being available
// to multiple jumping units. Range: 0.5–1.0
constexpr double GAMMA_VF = 0.6;

// V̂₁*: specific hole free volume required for a CO₂ molecule to jump
// For CO₂: ~0.589 cm³/g = 5.89e-4 m³/kg (from polymer literature)
// For heavy oil (PAG): ~0.750–0.900 cm³/g
constexpr double VSTAR_CO2 = 5.89e-4; // m³/kg
constexpr double VSTAR_OIL = 8.50e-4; // m³/kg

// Occupied volume per segment (from close-packed hard-sphere reference)
// Used as the "no-free-volume" reference density
constexpr double CP_FRAC = 0.74048;   // FCC close-packing fraction

// ============================================================================
// Viscosity model parameters  [PLACEHOLDERS]
// ============================================================================

// Free-volume viscosity model: μ = μ₀ × exp(γ·V̂*/V̂_FH)
// μ₀_i: pre-exponential viscosity for pure component i

// CO₂ viscosity at 0°C, low P: ~1.4e-5 Pa·s
// At high density (liquid-like): pre-factor for free-volume model
constexpr double MU0_CO2 = 1.2e-5;  // Pa·s

// SW68 (ISO VG 68) oil viscosity at 40°C: ~68 mPa·s = 6.8e-2 Pa·s
// Pre-factor estimated from free-volume back-calculation
constexpr double MU0_OIL = 4.5e-2;  // Pa·s

// Weighting exponent for mixture V̂* mixing rule
constexpr double VSTAR_MIX_EXP = 1.0;  // linear mixing

// ============================================================================
// Flow & Geometry parameters  [PLACEHOLDERS — replace with CFD input]
// ============================================================================

// L: characteristic length scale of the flow channel (e.g. tube diameter)
// Typical values: 0.001–0.1 m for heat exchanger channels
constexpr double L_CHAR = 0.01;   // m

// u: characteristic flow velocity (from CFD or bulk flow rate)
// Typical values: 0.1–10 m/s for refrigeration/heat-pump systems
constexpr double U_CHAR = 1.0;    // m/s

// ============================================================================
// Sherwood number correlation constants  [PLACEHOLDERS]
// ============================================================================

// Sh = C × Re^a × Sc^b
// Laminar boundary layer (flat plate): C=0.664, a=0.5, b=0.333
// Turbulent boundary layer (flat plate): C=0.037, a=0.8, b=0.333
// Internal laminar pipe flow: Sh ≈ 3.66 (constant)
// Internal turbulent pipe flow: C=0.023, a=0.8, b=0.333
//
// Default: laminar boundary layer (conservative estimate for oil-refrigerant)

constexpr double C_SH = 0.664;
constexpr double A_SH = 0.5;
constexpr double B_SH = 0.333;

// Minimum Sh to prevent division by zero in stagnant regions
constexpr double SH_MIN = 2.0;

// ============================================================================
// Unit conversions
// ============================================================================
constexpr double M3_PER_CM3  = 1.0e-6;   // cm³ → m³
constexpr double KG_PER_G    = 1.0e-3;   // g → kg

// ============================================================================
// T-dependent k_ij
// ============================================================================
double compute_kij(double T) {
    return A_k * T + B_k;
}

void update_eij_for_T(add_args& args, double T) {
    double kij = compute_kij(T);
    double si3  = args.s[0]*args.s[0]*args.s[0];
    double sj3  = args.s[1]*args.s[1]*args.s[1];
    double sij3 = args.s_ij[1]*args.s_ij[1]*args.s_ij[1];
    double eij  = (1.0 - kij) * std::sqrt(si3 * sj3) / sij3
                * std::sqrt(args.e[0] * args.e[1]);
    args.e_ij[1] = eij;
    args.e_ij[2] = eij;
}

// ============================================================================
// 1. Mixture density ρ_m (kg/m³) from SAFT EOS
// ============================================================================

/** Compute the mass density (kg/m³) of the CO2+oil liquid mixture at (T,P,x1).
 *  Uses SAFT molar density + molecular weight conversion.
 */
double compute_rho_mass(double T, double P, double x1, add_args& args) {
    std::vector<double> x = {x1, 1.0 - x1};
    double rho_molar = saft_rho_cpp(T, P, x, LIQUID, args);  // mol/m³

    // Convert to mass density: ρ_mass = ρ_molar × M̄
    double M_avg = x1 * MW_CO2 + (1.0 - x1) * MW_OIL;  // g/mol
    return rho_molar * M_avg * KG_PER_G;  // kg/m³
}

// ============================================================================
// 2. Free Volume Fraction  V̂_FH  (m³/kg)
// ============================================================================

/** Compute the average free volume fraction per unit mass.
 *
 *  Uses the SAFT packing fraction η:
 *    η = (π/6) × Σ(x_i × m_i × d_i³) × ρ_molar
 *    V_FH = (1 - η) / ρ_molar   (in m³/mol)
 *
 *  Converted to specific free volume (m³/kg):
 *    V̂_FH = V_FH / M̄
 */
double compute_VFH(double T, double rho_molar, double x1, add_args& args) {
    // Packing fraction from SAFT
    std::vector<double> x = {x1, 1.0 - x1};
    double eta = rho_to_eta(T, rho_molar, x, args);

    // Guard: η should be between 0 and ~0.74 (close-packed limit)
    if (eta >= CP_FRAC) eta = CP_FRAC * 0.999;
    if (eta <= 0.0)    eta = 0.001;

    // Free volume: 1 - η fraction of total volume
    double V_FH_molar = (1.0 - eta) / rho_molar;  // m³/mol
    double M_avg = x1 * MW_CO2 + (1.0 - x1) * MW_OIL;  // g/mol
    return V_FH_molar / (M_avg * KG_PER_G);  // m³/kg
}

// ============================================================================
// 3. Dynamic Diffusion Coefficient  D₁₂  (m²/s)
// ============================================================================

/** Compute the binary diffusion coefficient using free-volume theory.
 *
 *  Vrentas-Duda formulation:
 *    D₁₂ = D₀ × exp( -γ × V̂₁* / V̂_FH )
 *
 *  where V̂_FH is the average specific free volume of the mixture,
 *  and V̂₁* is the critical hole free volume for CO₂.
 */
double compute_D12(double T, double rho_molar, double x1, add_args& args) {
    double VFH = compute_VFH(T, rho_molar, x1, args);
    double exponent = -GAMMA_VF * VSTAR_CO2 / VFH;
    // Clamp the exponent to avoid numerical overflow
    if (exponent < -50.0) exponent = -50.0;
    if (exponent >  10.0) exponent =  10.0;
    return D0_CO2 * std::exp(exponent);
}

// ============================================================================
// 4. Mixture Viscosity  μ_m  (Pa·s)
// ============================================================================

/** Compute mixture viscosity using free-volume theory.
 *
 *  Consistent with the diffusion model: both use the same V̂_FH.
 *  Mixing rule for the critical free volume V̂*:
 *    V̂*_mix = w₁ × V̂*₁ + w₂ × V̂*₂
 *
 *  Free-volume viscosity model (Cohen-Turnbull / Doolittle):
 *    μ = μ₀ × exp( γ × V̂*_mix / V̂_FH )
 *
 *  The pure-component pre-factors μ₀ are blended using log-mixing:
 *    ln(μ₀_mix) = w₁ × ln(μ₀₁) + w₂ × ln(μ₀₂)
 */
double compute_viscosity(double T, double rho_molar, double x1, double w1,
                         add_args& args) {
    double VFH = compute_VFH(T, rho_molar, x1, args);
    if (VFH < 1e-20) VFH = 1e-20;

    double w2 = 1.0 - w1;

    // Mixing rule for V̂* and μ₀
    double Vstar_mix = w1 * VSTAR_CO2 + w2 * VSTAR_OIL;
    double ln_mu0_mix = w1 * std::log(MU0_CO2) + w2 * std::log(MU0_OIL);
    double mu0_mix = std::exp(ln_mu0_mix);

    double exponent = GAMMA_VF * Vstar_mix / VFH;
    if (exponent > 50.0) exponent = 50.0;
    return mu0_mix * std::exp(exponent);
}

// ============================================================================
// 5. Dimensionless Numbers and τ_ms
// ============================================================================

struct TransportResult {
    double rho_mass;   // kg/m³
    double D12;        // m²/s
    double mu_mix;     // Pa·s
    double Re;         // —
    double Sc;         // —
    double Sh;         // —
    double tau_ms;     // s
};

/** Compute the full transport chain for one (T,P,x1,w1) point.
 *
 *  The mass fraction w1 is needed for viscosity mixing rules.
 *  x1 and w1 are related by MW mapping, but both are provided for convenience.
 */
TransportResult compute_transport(double T, double P, double x1, double w1,
                                   add_args& args) {
    TransportResult tr{};

    // Step 1: Mass density
    tr.rho_mass = compute_rho_mass(T, P, x1, args);
    std::vector<double> x = {x1, 1.0 - x1};
    double rho_molar = saft_rho_cpp(T, P, x, LIQUID, args);

    // Step 2: Diffusion coefficient
    tr.D12 = compute_D12(T, rho_molar, x1, args);

    // Step 3: Viscosity
    tr.mu_mix = compute_viscosity(T, rho_molar, x1, w1, args);

    // Step 4: Dimensionless numbers
    tr.Re = tr.rho_mass * U_CHAR * L_CHAR / tr.mu_mix;
    tr.Sc = tr.mu_mix / (tr.rho_mass * tr.D12);
    tr.Sh = C_SH * std::pow(tr.Re, A_SH) * std::pow(tr.Sc, B_SH);
    if (tr.Sh < SH_MIN) tr.Sh = SH_MIN;

    // Step 5: Characteristic mass transfer time
    tr.tau_ms = L_CHAR * L_CHAR / (tr.D12 * tr.Sh);

    return tr;
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "=== Step 2: CO2 + SW68 Mesoscopic Transport (tau_ms) ===" << std::endl;
    std::cout << "\n--- Physical Parameters ---" << std::endl;
    std::cout << "D0_CO2   = " << D0_CO2   << " m^2/s" << std::endl;
    std::cout << "GAMMA_VF  = " << GAMMA_VF  << std::endl;
    std::cout << "VSTAR_CO2 = " << VSTAR_CO2 << " m^3/kg" << std::endl;
    std::cout << "VSTAR_OIL = " << VSTAR_OIL << " m^3/kg" << std::endl;
    std::cout << "MU0_CO2  = " << MU0_CO2  << " Pa.s" << std::endl;
    std::cout << "MU0_OIL  = " << MU0_OIL  << " Pa.s" << std::endl;
    std::cout << "\n--- Flow Parameters ---" << std::endl;
    std::cout << "L_char   = " << L_CHAR   << " m" << std::endl;
    std::cout << "u_char   = " << U_CHAR   << " m/s" << std::endl;
    std::cout << "C_sh     = " << C_SH  << " (laminar flat plate)" << std::endl;
    std::cout << "a_sh     = " << A_SH  << std::endl;
    std::cout << "b_sh     = " << B_SH  << std::endl;

    // Setup SAFT
    add_args args;
    args.m = {m_CO2, m_OIL};  args.s = {s_CO2, s_OIL};  args.e = {e_CO2, e_OIL};
    args.lambda_a = {la_CO2, la_OIL};  args.lambda_r = {lr_CO2, lr_OIL};
    init_saft_cross_params(args);

    // Grid (same as Step 1)
    constexpr double T_min = 283.15, T_max = 343.15, dT = 1.0;
    constexpr double P_min = 3.0e6,  P_max = 15.0e6,  dP = 0.1e6;
    constexpr int nT = static_cast<int>((T_max - T_min) / dT + 0.5) + 1;
    constexpr int nP = static_cast<int>((P_max - P_min) / dP + 0.5) + 1;
    std::cout << "\nGrid: " << nT << " x " << nP << " = " << nT*nP << " nodes" << std::endl;

    // Read w1 LUT from Step 1
    std::ifstream fin("weq_lut.csv");
    if (!fin.is_open()) {
        std::cerr << "ERROR: weq_lut.csv not found. Run Step 1 (calc_weq_lut) first." << std::endl;
        return 1;
    }
    std::string header;
    std::getline(fin, header);

    // We only need T, P, x1, w1 from the Step 1 output
    // The file has columns: T_K,P_Pa,x1,w1_eq,rho_L_molm3,phi1_V,phi1_L,iter,converged
    std::vector<double> x1_arr, w1_arr;
    x1_arr.reserve(nT * nP);
    w1_arr.reserve(nT * nP);
    std::string line;
    while (std::getline(fin, line)) {
        // Parse x1 (3rd field) and w1 (4th field)
        size_t p1 = line.find(',');
        size_t p2 = line.find(',', p1 + 1);
        size_t p3 = line.find(',', p2 + 1);
        size_t p4 = line.find(',', p3 + 1);
        if (p4 == std::string::npos) continue;
        double val_w1 = std::stod(line.substr(p3 + 1, p4 - p3 - 1));
        w1_arr.push_back(val_w1);
        // x1 is between p2 and p3
        double val_x1 = std::stod(line.substr(p2 + 1, p3 - p2 - 1));
        x1_arr.push_back(val_x1);
    }
    fin.close();

    if ((int)x1_arr.size() != nT * nP) {
        std::cerr << "ERROR: Expected " << nT*nP << " rows, got " << x1_arr.size() << std::endl;
        return 1;
    }
    std::cout << "Loaded " << x1_arr.size() << " points from weq_lut.csv" << std::endl;

    // Output
    std::ofstream fout("taums_lut.csv");
    fout << std::scientific << std::setprecision(10);
    fout << "T_K,P_Pa,x1,w1_eq,rho_kgm3,D12_m2s,mu_Pas,Re,Sc,Sh,tau_ms_s" << std::endl;

    int n_ok = 0, n_skip = 0;
    TransportResult tr_min{0,1e10,0,0,0,0,1e10};
    TransportResult tr_max{0,0,0,0,0,0,0};

    for (int iT = 0; iT < nT; iT++) {
        double T = T_min + iT * dT;
        update_eij_for_T(args, T);

        for (int iP = 0; iP < nP; iP++) {
            double P = P_min + iP * dP;
            int idx = iT * nP + iP;
            double x1 = x1_arr[idx];
            double w1 = w1_arr[idx];

            // Skip full-miscibility points (x1 ≈ 1 gives unrealistic transport)
            if (x1 > 0.999) {
                fout << T << "," << P << "," << x1 << "," << w1 << ","
                     << "nan,nan,nan,nan,nan,nan,nan" << std::endl;
                n_skip++;
                continue;
            }

            try {
                TransportResult tr = compute_transport(T, P, x1, w1, args);

                fout << T << "," << P << "," << x1 << "," << w1 << ","
                     << tr.rho_mass << "," << tr.D12 << "," << tr.mu_mix << ","
                     << tr.Re << "," << tr.Sc << "," << tr.Sh << ","
                     << tr.tau_ms << std::endl;

                // Track min/max for reporting
                if (tr.tau_ms < tr_min.tau_ms) tr_min = tr;
                if (tr.tau_ms > tr_max.tau_ms) tr_max = tr;
                n_ok++;
            } catch (const std::exception&) {
                fout << T << "," << P << "," << x1 << "," << w1 << ","
                     << "nan,nan,nan,nan,nan,nan,nan" << std::endl;
                n_skip++;
            }
        }
        if ((iT + 1) % 10 == 0) {
            std::cout << "  T = " << T << " K  (" << iT+1 << "/" << nT
                      << ")  ok=" << n_ok << "  skip=" << n_skip << std::endl;
        }
    }

    fout.close();

    std::cout << "\n==============================================" << std::endl;
    std::cout << "  Computed:    " << n_ok  << std::endl;
    std::cout << "  Skipped:     " << n_skip << std::endl;
    std::cout << "  tau_ms range:" << std::endl;
    if (n_ok > 0) {
        std::cout << "    min: " << tr_min.tau_ms << " s"
                  << "  (T=" << tr_min.rho_mass << ")" << std::endl;
        // Reconstruct T,P for min point
        std::cout << "    max: " << tr_max.tau_ms << " s" << std::endl;
    }
    std::cout << "  Output: taums_lut.csv" << std::endl;
    std::cout << "==============================================" << std::endl;

    return 0;
}
