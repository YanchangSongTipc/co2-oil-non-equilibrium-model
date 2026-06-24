---
name: CO2-Non-Equilibrium-Project
description: SAFT-VR-Mie CO2-oil dissolution model — code layout, parameters, compute chain, grid specs
type: project
---

# CO2 + SW68 Lubricant Oil Non-Equilibrium Dissolution Model

## Project Goal
Validate a non-equilibrium dissolution model for CO2 in polyol-ester (SW68) lubricating oil,
using SAFT-VR-Mie equation of state. The output feeds a macroscopic 1-D flow/heat-transfer model.

## Code Layout

| File | Purpose |
|------|---------|
| `src/saft_mie.h` / `src/saft_mie.cpp` | SAFT-VR-Mie EOS core (from Reaktoro project), with autodiff |
| `lib/autodiff/` | Header-only forward-mode automatic differentiation library |
| `src/calc_weq_lut.cpp` | **Step 1**: equilibrium solubility LUT generator |
| `src/calc_taums.cpp` | **Step 2**: mesoscopic transport tau_ms LUT generator |
| `src/test.cpp` | Pure CO2 verification test |
| `CMakeLists.txt` | Build system |
| `scripts/plot_weq.py` | Solubility visualisation |
| `scripts/plot_taums.py` | Transport visualisation + sensitivity analysis |
| `docs/what_needed.md` | Original task specification from the research report |

## Component SAFT Parameters

### CO2 (component 0)
```
MW = 44.0
m  = 1.2625
sigma = 3.5183 Angstrom
e/k  = 261.62 K
lambda_r = 21.245
lambda_a = 5.687
```

### SW68 Lubricant Oil (component 1)
```
MW = 634.0
m  = 12.087
sigma = 4.391 Angstrom
e/k  = 405.77 K
lambda_r = 18.0
lambda_a = 6.0
```

### Temperature-dependent k_ij
```
k_ij(T) = A * T + B
A = -6.719e-7 * 634 + 5.339e-4 = 1.078e-4
B =  2.979e-4 * 634 - 0.217    = -0.02813
```

k_ij(283.15 K) = 0.0024
k_ij(343.15 K) = 0.0089

## Compute Grid
- Temperature: 283.15 – 343.15 K, dT = 1 K (61 points)
- Pressure: 3.0 – 15.0 MPa, dP = 0.1 MPa (121 points)
- Total: 7,381 nodes per LUT

## Build Commands (Ubuntu 22.04)
```bash
sudo apt install build-essential g++ python3-pip
pip install numpy pandas matplotlib scipy

# Compile (autodiff headers in lib/, sources in src/)
g++ -std=c++17 -O2 -I lib src/calc_weq_lut.cpp src/saft_mie.cpp -o calc_weq_lut
g++ -std=c++17 -O2 -I lib src/calc_taums.cpp  src/saft_mie.cpp -o calc_taums

# Run (order matters)
./calc_weq_lut          # produces weq_lut.csv
./calc_taums            # reads weq_lut.csv, produces taums_lut.csv
python scripts/plot_weq.py      # solubility figures
python scripts/plot_taums.py    # transport figures
```

## Key Design Decisions
1. **Solver**: damped secant method with bracket detection. Pure secant oscillates near P_sat;
   required step-size damping (max 30-40% of current x1) for stability.
2. **Phase boundary handling**: Once an isotherm enters full-miscibility (x1 > 0.999),
   all higher P at that T are locked to full-miscibility. Prevents oscillation between
   two-phase and full-misc branches.
3. **Free volume theory**: V̂_FH computed from SAFT packing fraction eta via `rho_to_eta()`.
4. **Viscosity model**: Free-volume theory consistent with diffusion model (shared V̂_FH).
5. **All placeholder parameters** (D0, gamma, Vstar_*, mu0_*, L, u, C,a,b) are in
   `calc_taums.cpp` lines 42-80. Replace with experimental values when available.

## Current Status
- Step 1: 7,381/7,381 converged (100%), w1 range 0.082-0.999
- Step 2: 4,718/7,381 computed (rest are full-miscibility), tau_ms range 16-32 s
- tau_ms values are ~1000x too high due to placeholder D0 and L; will calibrate later

## Fixed Issues (do not reintroduce)
1. DO NOT uncomment `printf("g_HS: ...")` at saft_mie.cpp line 209 — floods stdout, kills performance
2. DO NOT use pure secant without damping — oscillates at P_sat
3. DO NOT remove the full-miscibility lock (`entered_full_misc` flag) — causes alternating oscillations
4. The `\\r` in `\\rho` in matplotlib titles causes LaTeX parse errors on some systems — use raw strings `r'$\\rho$'`
