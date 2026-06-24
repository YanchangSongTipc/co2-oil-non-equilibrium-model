#!/usr/bin/env python3
"""
Visualization of CO2 + SW68 equilibrium solubility LUT.
Standard phase-equilibrium convention:
  - X-axis: mass fraction w1_eq (solubility)
  - Y-axis: pressure P (MPa)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import os

# ── Load data ────────────────────────────────────────────────────────────────
# Run from project root regardless of where this script is located
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

df = pd.read_csv('weq_lut.csv')
print(f"Loaded {len(df)} rows, {df['converged'].sum()} converged")

# Pivot into 2D grid
T_arr = np.sort(df['T_K'].unique())
P_arr = np.sort(df['P_Pa'].unique()) / 1e6  # MPa
nT, nP = len(T_arr), len(P_arr)
W = df['w1_eq'].values.reshape(nT, nP)
X1 = df['x1'].values.reshape(nT, nP)
RhoL = df['rho_L_molm3'].values.reshape(nT, nP)
TT, PP = np.meshgrid(T_arr, P_arr, indexing='ij')

print(f"Grid: {nT} T x {nP} P")
print(f"T range: {T_arr[0]:.1f} -- {T_arr[-1]:.1f} K")
print(f"P range: {P_arr[0]:.1f} -- {P_arr[-1]:.1f} MPa")
print(f"w1_eq range: {W.min():.4f} -- {W.max():.4f}")
print(f"x1 range:    {X1.min():.4f} -- {X1.max():.4f}")

# ── Global style ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'legend.fontsize': 9, 'figure.dpi': 150,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

# Colour palette for isotherms
iso_T_all = [283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 323.15, 333.15, 343.15]
iso_T_main = [283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15]
colors_main = plt.cm.plasma(np.linspace(0, 0.95, len(iso_T_main)))
colors_all = plt.cm.plasma(np.linspace(0, 0.95, len(iso_T_all)))

# ── 1. P-w Isotherms (primary figure) ────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: mass fraction
for i, T_sel in enumerate(iso_T_all):
    iT = np.argmin(np.abs(T_arr - T_sel))
    axes[0].plot(W[iT, :], P_arr, color=colors_all[i], lw=1.5,
                 label=f'{T_sel:.0f} K' if T_sel in iso_T_main else None)
axes[0].set_xlabel(r'CO$_2$ Mass Fraction  $w_{1,eq}$  (kg/kg)')
axes[0].set_ylabel('Pressure  $P$  (MPa)')
axes[0].set_title('(a)  Solubility Isotherms  $P$–$w_{1,eq}$', fontsize=13)
axes[0].legend(loc='lower right', ncol=2, fontsize=9)
axes[0].grid(True, alpha=0.3)
axes[0].set_xlim(left=0)

# Right: mole fraction
for i, T_sel in enumerate(iso_T_all):
    iT = np.argmin(np.abs(T_arr - T_sel))
    axes[1].plot(X1[iT, :], P_arr, color=colors_all[i], lw=1.5,
                 label=f'{T_sel:.0f} K' if T_sel in iso_T_main else None)
axes[1].set_xlabel(r'CO$_2$ Mole Fraction  $x_1$  (mol/mol)')
axes[1].set_ylabel('Pressure  $P$  (MPa)')
axes[1].set_title('(b)  Solubility Isotherms  $P$–$x_1$', fontsize=13)
axes[1].legend(loc='lower right', ncol=2, fontsize=9)
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim(left=0)

fig.suptitle('CO$_2$ (1) + SW68 (2):  Isothermal Solubility Curves\n'
             r'$k_{ij}(T) = A \cdot T + B$,  $A=1.078\times10^{-4}$,  $B=-0.0281$',
             fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('fig_isotherms.png', dpi=200)
plt.close()
print("[OK] fig_isotherms.png")

# ── 2. T-P Contour Map ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
cf = ax.contourf(TT, PP, W, levels=30, cmap='RdYlBu_r')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label(r'$w_{1,eq}$  (kg CO$_2$ / kg mixture)')
ct = ax.contour(TT, PP, W, levels=[0.1, 0.2, 0.3, 0.5, 0.7, 0.9],
                colors='black', linewidths=0.6, linestyles='--')
ax.clabel(ct, fmt='%.1f', fontsize=8)
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  $w_{1,eq}(T,P)$  Filled Contour', fontsize=14)
fig.savefig('fig_contour.png', dpi=200)
plt.close()
print("[OK] fig_contour.png")

# ── 3. 3D Surface ────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(TT, PP, W, cmap=cm.viridis, linewidth=0, antialiased=True, alpha=0.92)
ax.set_xlabel('T (K)')
ax.set_ylabel('P (MPa)')
ax.set_zlabel(r'$w_{1,eq}$')
ax.set_title(r'CO$_2$ + SW68  $w_{1,eq}(T,P)$  3D Surface', fontsize=14, pad=15)
cbar = fig.colorbar(surf, ax=ax, shrink=0.55, pad=0.08)
cbar.set_label(r'$w_{1,eq}$')
ax.view_init(elev=22, azim=-55)
fig.savefig('fig_3d_surface.png', dpi=200)
plt.close()
print("[OK] fig_3d_surface.png")

# ── 4. P-w1 Contour (phase-diagram style) ─────────────────────────────────────
# Build P-w grid for each T and overlay
fig, ax = plt.subplots(figsize=(9, 8))
for i, T_sel in enumerate(iso_T_all):
    iT = np.argmin(np.abs(T_arr - T_sel))
    w = W[iT, :]
    # Filter out jumps to w1≈1 (full miscibility) for smoother curves
    # Keep points where w1 increases monotonically, stop at first jump
    mask = np.ones(len(w), dtype=bool)
    for j in range(1, len(w)):
        if w[j] - w[j-1] > 0.3:  # jump > 0.3 → phase transition
            mask[j:] = False
            break
    w_plot = w[mask]
    p_plot = P_arr[mask]
    ax.plot(w_plot, p_plot, color=colors_all[i], lw=1.5,
            label=f'{T_sel:.0f} K' if T_sel in iso_T_main else None,
            alpha=0.85)

ax.set_xlabel(r'CO$_2$ Mass Fraction  $w_{1,eq}$  (kg/kg)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  $P$–$w_{1,eq}$  Phase Diagram', fontsize=14)
ax.legend(loc='lower right', ncol=2, fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig('fig_phase_diagram.png', dpi=200)
plt.close()
print("[OK] fig_phase_diagram.png")

# ── 5. Phase Regime Map ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
regime = np.where(X1 >= 0.999, 2, 1)
cf = ax.pcolormesh(TT, PP, regime, cmap=plt.cm.RdYlBu, shading='auto',
                   vmin=0.5, vmax=2.5)
cbar = fig.colorbar(cf, ax=ax, ticks=[1, 2])
cbar.ax.set_yticklabels([r'Two-Phase Equilibrium  ($x_1 < 0.999$)',
                          r'Full Miscibility  ($x_1 \approx 1$)'])
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  Phase Regime Map  ($T$–$P$)', fontsize=14)
fig.savefig('fig_regime.png', dpi=200)
plt.close()
print("[OK] fig_regime.png")

# ── 6. Density Map ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
cf = ax.contourf(TT, PP, RhoL / 1000, levels=30, cmap='YlOrRd')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label(r'Liquid Density  $\rho_L$  ($\times 10^3$ mol/m$^3$)')
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  Liquid Mixture Density  $\rho_L(T,P)$', fontsize=14)
fig.savefig('fig_density.png', dpi=200)
plt.close()
print("[OK] fig_density.png")

# ── 7. Combined Report Figure ────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 13))
gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.35)

# (a) P-w Isotherms
ax1 = fig.add_subplot(gs[0, 0])
for i, T_sel in enumerate(iso_T_all):
    iT = np.argmin(np.abs(T_arr - T_sel))
    w = W[iT, :]
    mask = np.ones(len(w), dtype=bool)
    for j in range(1, len(w)):
        if w[j] - w[j-1] > 0.3:
            mask[j:] = False; break
    ax1.plot(w[mask], P_arr[mask], color=colors_all[i], lw=1.2,
             label=f'{T_sel:.0f} K' if T_sel in iso_T_main else None)
ax1.set_xlabel(r'$w_{1,eq}$ (kg/kg)')
ax1.set_ylabel(r'$P$ (MPa)')
ax1.set_title(r'(a)  $P$–$w_{1,eq}$ Isotherms')
ax1.legend(fontsize=7, ncol=2, loc='lower right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(left=0)

# (b) T-P Contour
ax2 = fig.add_subplot(gs[0, 1:])
cf2 = ax2.contourf(TT, PP, W, levels=30, cmap='RdYlBu_r')
cbar2 = fig.colorbar(cf2, ax=ax2, shrink=0.8)
cbar2.set_label(r'$w_{1,eq}$', fontsize=9)
ct2 = ax2.contour(TT, PP, W, levels=[0.1, 0.2, 0.3, 0.5, 0.7, 0.9],
                  colors='black', linewidths=0.5, linestyles='--')
ax2.clabel(ct2, fmt='%.1f', fontsize=7)
ax2.set_xlabel(r'$T$ (K)')
ax2.set_ylabel(r'$P$ (MPa)')
ax2.set_title(r'(b)  $w_{1,eq}(T,P)$  Contour Map')

# (c) Regime Map
ax3 = fig.add_subplot(gs[1, 0])
cf3 = ax3.pcolormesh(TT, PP, regime, cmap=plt.cm.RdYlBu, shading='auto',
                     vmin=0.5, vmax=2.5)
cbar3 = fig.colorbar(cf3, ax=ax3, ticks=[1, 2], shrink=0.7)
cbar3.ax.set_yticklabels(['Two-Phase', 'Full-Misc.'], fontsize=8)
ax3.set_xlabel(r'$T$ (K)')
ax3.set_ylabel(r'$P$ (MPa)')
ax3.set_title(r'(c)  Phase Regime')

# (d) P-x1 Isotherms
ax4 = fig.add_subplot(gs[1, 1])
for i, T_sel in enumerate(iso_T_main):
    iT = np.argmin(np.abs(T_arr - T_sel))
    ax4.plot(X1[iT, :], P_arr, color=colors_main[i], lw=1.5, label=f'{T_sel:.0f} K')
ax4.set_xlabel(r'$x_1$ (mol/mol)')
ax4.set_ylabel(r'$P$ (MPa)')
ax4.set_title(r'(d)  $P$–$x_1$ Isotherms')
ax4.legend(fontsize=8, ncol=2, loc='lower right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(left=0)

# (e) Density
ax5 = fig.add_subplot(gs[1, 2])
cf5 = ax5.contourf(TT, PP, RhoL / 1000, levels=30, cmap='YlOrRd')
cbar5 = fig.colorbar(cf5, ax=ax5, shrink=0.7)
cbar5.set_label(r'$\rho_L$ ($\times 10^3$ mol/m$^3$)', fontsize=9)
ax5.set_xlabel(r'$T$ (K)')
ax5.set_ylabel(r'$P$ (MPa)')
ax5.set_title(r'(e)  Density $\rho_L(T,P)$')

fig.suptitle(r'CO$_2$ (1) + SW68 Lubricant Oil (2) — SAFT-VR-Mie Phase Equilibrium LUT' '\n'
             r'$k_{ij}(T) = A\cdot T + B$,  $A = 1.078\times 10^{-4}$,  $B = -0.0281$',
             fontsize=15, fontweight='bold', y=1.01)
fig.savefig('fig_report.png', dpi=200)
plt.close()
print("[OK] fig_report.png (combined 5-panel figure)")

# ── 8. Data Summary Table ────────────────────────────────────────────────────
sample_T = [283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15]
sample_P = [3.0, 5.0, 8.0, 12.0, 15.0]

print("\n" + "=" * 95)
print("  Key Data:  w1_eq  (CO2 mass fraction in liquid)")
print("=" * 95)
header = f"{'T (K)':>8}"
for p in sample_P:
    header += f"  P={p:.0f} MPa"
print(header)
print("-" * 95)
for T_sel in sample_T:
    iT = np.argmin(np.abs(T_arr - T_sel))
    row = f"{T_sel:8.1f}"
    for P_sel in sample_P:
        iP = np.argmin(np.abs(P_arr - P_sel))
        row += f"  {W[iT, iP]:11.6f}"
    print(row)
print("-" * 95)
print("  Note: 0.998561 = full-miscibility regime (x1 >= 0.999)")
print(f"         CO2 P_sat(283 K) ~ 4.5 MPa;  Tc(CO2) = 304.1 K")
print("=" * 95)

# ── Summary ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("All figures saved:")
files = ['fig_isotherms.png', 'fig_phase_diagram.png', 'fig_contour.png',
         'fig_3d_surface.png', 'fig_regime.png', 'fig_density.png', 'fig_report.png']
for f in files:
    try:
        size_kb = os.path.getsize(f) / 1024
        print(f"  {f:30s}  ({size_kb:.0f} KB)")
    except:
        print(f"  {f:30s}  (MISSING)")
print("=" * 60)
