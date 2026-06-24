---
name: Implementation Feedback
description: Solver design lessons, numerical gotchas, and validation notes from this implementation
type: feedback
---

# Implementation Feedback

## Rule 1: Damped secant, not pure secant
**Why**: The phase equilibrium function f(x1) = phi_V - x1*phi_L(x1) has two roots
near P_sat (physical x1<1 and trivial x1=1). Pure secant oscillates between them.
**How to apply**: Always use step-size damping (max_step = 0.3*(x_sec0 + x_sec1)) and
bracket detection. The `src/calc_weq_lut.cpp` solver is the reference implementation.

## Rule 2: Full-miscibility lock-in for isotherms
**Why**: Once P exceeds the CO2 saturation pressure at a given T, the physical system
enters single-phase (full miscibility). Without locking, the solver oscillates between
finding two-phase solutions (metastable) and full-miscibility at alternating P steps.
**How to apply**: Use the `entered_full_misc` boolean flag — once x1 > 0.999, skip
all higher P at that T. Implemented in `calc_weq_lut.cpp` main loop.

## Rule 3: Debug print in SAFT must stay commented
**Why**: `printf("g_HS: ...")` at saft_mie.cpp line 209 prints at every radial distribution
function evaluation (hot path). This floods output with millions of lines and kills performance.
**How to apply**: Never uncomment line 209. If debugging SAFT, add targeted prints with
conditional guards.

## Rule 4: Python `\r` in matplotlib LaTeX
**Why**: On Windows GBK locale, `\r` in `\rho` is interpreted as carriage return.
On Ubuntu with UTF-8 locale, this may or may not cause issues depending on matplotlib backend.
**How to apply**: Always use raw strings: `r'$\rho_L(T,P)$'` instead of `'$\rho_L(T,P)$'`.

## Rule 5: Proxy configuration for Git
**Why**: Ubuntu machine may be behind institutional firewall/proxy. Git operations
fail without proxy configuration.
**How to apply**: Set `http.proxy` and `https.proxy` in git config, or use
`export HTTPS_PROXY=http://proxy:port` before gh commands.

## Rule 6: CSV reading assumes exact grid dimensions
**Why**: calc_taums.cpp reads weq_lut.csv and expects exactly nT*nP = 7381 rows.
If the grid changes, this will fail with a clear error message.
**How to apply**: Always run Step 1 before Step 2. If changing the grid, update
both source files consistently.

## Rule 7: Cross-platform: no Windows-specific code
**Why**: User is migrating from Windows 11 to Ubuntu 22.04. Code must compile and
run on both.
**How to apply**: Use standard C++17, forward-slash paths in Python, no `#pragma`
that is MSVC-only. Current codebase is clean — only `#pragma GCC optimize(3)` in
test.cpp which is GCC-compatible on Linux.
