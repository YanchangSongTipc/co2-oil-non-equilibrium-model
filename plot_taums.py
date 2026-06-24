#!/usr/bin/env python3
"""
Visualization of tau_ms (mass transfer relaxation time) LUT.
Plots tau_ms as a function of T and P for the CO2+SW68 system.
Also includes sensitivity analysis for the key unknown parameters.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Load data
df = pd.read_csv('taums_lut.csv')
T_arr = np.sort(df['T_K'].unique())
P_arr = np.sort(df['P_Pa'].unique()) / 1e6
nT, nP = len(T_arr), len(P_arr)

# Reshape into grids — use NaN for skipped (full-misc) points
def reshape_col(col):
    data = np.full(nT * nP, np.nan)
    raw = df[col].values
    for i in range(len(raw)):
        data[i] = raw[i]
    return data.reshape(nT, nP)

Tau = reshape_col('tau_ms_s') * 1000  # convert to ms
D12v = reshape_col('D12_m2s')
Mu   = reshape_col('mu_Pas')
TT, PP = np.meshgrid(T_arr, P_arr, indexing='ij')

ok_mask = ~np.isnan(Tau)
print(f"Valid points: {np.sum(ok_mask)} / {nT*nP}")
if np.sum(ok_mask) > 0:
    print(f"tau_ms: {np.nanmin(Tau):.1f} – {np.nanmax(Tau):.1f} ms")
    print(f"D12:    {np.nanmin(D12v):.2e} – {np.nanmax(D12v):.2e} m2/s")
    print(f"mu:     {np.nanmin(Mu):.4f} – {np.nanmax(Mu):.4f} Pa.s")

plt.rcParams.update({
    'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 12,
    'figure.dpi': 150, 'savefig.bbox': 'tight',
})

# ── Colour palettes ──────────────────────────────────────────────────────────
iso_T = [283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15]
colors = plt.cm.plasma(np.linspace(0, 0.95, len(iso_T)))

# ── 1. tau_ms Contour Map ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
# Mask NaN regions
Tau_masked = np.ma.masked_invalid(Tau)
cf = ax.contourf(TT, PP, Tau_masked, levels=30, cmap='YlOrRd_r')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label(r'$\tau_{ms}$  (ms)')
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  Mass Transfer Relaxation Time  $\tau_{ms}(T,P)$', fontsize=14)
fig.savefig('fig_taums_contour.png', dpi=200)
plt.close()
print("[OK] fig_taums_contour.png")

# ── 2. D12 Contour Map ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
D12_masked = np.ma.masked_invalid(D12v)
cf = ax.contourf(TT, PP, D12_masked * 1e9, levels=30, cmap='viridis')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label(r'$D_{12}$  ($\times 10^{-9}$ m$^2$/s)')
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  Diffusion Coefficient  $D_{12}(T,P)$', fontsize=14)
fig.savefig('fig_d12_contour.png', dpi=200)
plt.close()
print("[OK] fig_d12_contour.png")

# ── 3. Viscosity Contour Map ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
Mu_masked = np.ma.masked_invalid(Mu)
cf = ax.contourf(TT, PP, Mu_masked * 1000, levels=30, cmap='YlOrRd_r')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label(r'$\mu_m$  (mPa$\cdot$s)')
ax.set_xlabel('Temperature  $T$  (K)')
ax.set_ylabel('Pressure  $P$  (MPa)')
ax.set_title(r'CO$_2$ + SW68:  Mixture Viscosity  $\mu_m(T,P)$', fontsize=14)
fig.savefig('fig_viscosity_contour.png', dpi=200)
plt.close()
print("[OK] fig_viscosity_contour.png")

# ── 4. Isothermal transport curves ───────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for i, T_sel in enumerate(iso_T):
    iT = np.argmin(np.abs(T_arr - T_sel))
    tau_row = Tau[iT, :]
    d12_row = D12v[iT, :]
    mu_row  = Mu[iT, :]
    # Filter NaN
    mask = ~np.isnan(tau_row)
    if np.sum(mask) == 0: continue
    label = f'{T_sel:.0f} K'
    axes[0,0].plot(P_arr[mask], tau_row[mask], color=colors[i], lw=1.5, label=label)
    axes[0,1].plot(P_arr[mask], d12_row[mask]*1e9, color=colors[i], lw=1.5, label=label)
    axes[1,0].plot(P_arr[mask], mu_row[mask]*1000, color=colors[i], lw=1.5, label=label)
    mask_w1 = mask & (np.abs(np.diff(np.concatenate([[0], tau_row]))) < 100)
    if np.sum(mask) > 0:
        w1_idx = np.argmin(np.abs(df[df['T_K']==T_arr[iT]]['w1_eq'].values[:nP] - 0.3))

for ax, ylab, title in [
    (axes[0,0], r'$\tau_{ms}$ (ms)', r'Mass Transfer Relaxation Time'),
    (axes[0,1], r'$D_{12}$ ($\times 10^{-9}$ m$^2$/s)', r'Diffusion Coefficient'),
    (axes[1,0], r'$\mu_m$ (mPa$\cdot$s)', r'Mixture Viscosity'),
]:
    ax.set_xlabel('Pressure  $P$  (MPa)')
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=8, ncol=2, loc='best')
    ax.grid(True, alpha=0.3)

# Bottom-right: tau_ms vs w1
ax = axes[1,1]
for i, T_sel in enumerate(iso_T):
    iT = np.argmin(np.abs(T_arr - T_sel))
    tau_row = Tau[iT, :]
    w1_row = df[df['T_K']==T_arr[iT]]['w1_eq'].values[:nP]
    mask = ~np.isnan(tau_row) & (w1_row < 0.999)
    if np.sum(mask) < 2: continue
    ax.plot(w1_row[mask], tau_row[mask], 'o-', color=colors[i], ms=3, lw=0.8,
            label=f'{T_sel:.0f} K')
ax.set_xlabel(r'$w_{1,eq}$  (kg/kg)')
ax.set_ylabel(r'$\tau_{ms}$  (ms)')
ax.set_title(r'$\tau_{ms}$ vs Solubility', fontsize=12)
ax.legend(fontsize=8, ncol=2, loc='best')
ax.grid(True, alpha=0.3)

fig.suptitle('CO$_2$ + SW68:  Mesoscopic Transport Properties\n'
             r'$D_0=5\times10^{-8}$ m$^2$/s,  $\gamma=0.6$,  $L=0.01$ m,  $u=1$ m/s  [PLACEHOLDER]',
             fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('fig_transport_panels.png', dpi=200)
plt.close()
print("[OK] fig_transport_panels.png")

# ── 5. Sensitivity analysis ──────────────────────────────────────────────────
# Show how tau_ms scales with D0 and L (the two most uncertain parameters)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Fix T=323K, P=8MPa as a representative point
iT_ref = np.argmin(np.abs(T_arr - 323.15))
iP_ref = np.argmin(np.abs(P_arr - 8.0))
tau_ref = Tau[iT_ref, iP_ref]
d12_ref = D12v[iT_ref, iP_ref] if not np.isnan(D12v[iT_ref, iP_ref]) else 2.8e-8

# Constants (matching C++ calc_taums.cpp)
D0_CO2_v = 5.0e-8  # m2/s

# Sensitivity to D0
D0_vals = np.logspace(-9, -5, 50)  # 1e-9 to 1e-5
tau_D0 = tau_ref * (5e-8 / D0_vals)  # scale relative to default D0=5e-8
axes[0].loglog(D0_vals, tau_D0, 'b-', lw=2)
axes[0].axvline(5e-8, color='red', ls='--', alpha=0.5, label='Current D0')
axes[0].axhline(1.0, color='green', ls='--', alpha=0.5, label='1 ms target')
axes[0].axhline(10.0, color='orange', ls='--', alpha=0.5, label='10 ms target')
axes[0].set_xlabel(r'$D_0$  (m$^2$/s)')
axes[0].set_ylabel(r'$\tau_{ms}$  (ms)')
axes[0].set_title(r'Sensitivity to $D_0$  (at T=323K, P=8MPa)')
axes[0].legend()
axes[0].grid(True, alpha=0.3, which='both')

# Sensitivity to L
L_vals = np.logspace(-4, -1, 50)  # 0.1 mm to 10 cm
tau_L = tau_ref * (L_vals / 0.01)**2  # tau ∝ L²
axes[1].loglog(L_vals * 1000, tau_L, 'b-', lw=2)  # x-axis in mm
axes[1].axvline(10, color='red', ls='--', alpha=0.5, label='Current L=10 mm')
axes[1].axhline(1.0, color='green', ls='--', alpha=0.5, label='1 ms target')
axes[1].axhline(100.0, color='orange', ls='--', alpha=0.5, label='100 ms target')
axes[1].set_xlabel(r'$L$  (mm)')
axes[1].set_ylabel(r'$\tau_{ms}$  (ms)')
axes[1].set_title(r'Sensitivity to $L$  (at T=323K, P=8MPa)')
axes[1].legend()
axes[1].grid(True, alpha=0.3, which='both')

fig.suptitle(r'CO$_2$ + SW68:  $\tau_{ms}$ Parameter Sensitivity',
             fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig('fig_sensitivity.png', dpi=200)
plt.close()
print("[OK] fig_sensitivity.png")

# ── 6. Combined report figure ────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# (a) tau_ms contour
ax1 = fig.add_subplot(gs[0, :2])
cf1 = ax1.contourf(TT, PP, Tau_masked, levels=30, cmap='YlOrRd_r')
cbar1 = fig.colorbar(cf1, ax=ax1, shrink=0.8)
cbar1.set_label(r'$\tau_{ms}$ (ms)', fontsize=9)
ax1.set_xlabel(r'$T$ (K)')
ax1.set_ylabel(r'$P$ (MPa)')
ax1.set_title(r'(a)  $\tau_{ms}(T,P)$  Mass Transfer Relaxation Time')

# (b) D12 contour
ax2 = fig.add_subplot(gs[0, 2])
cf2 = ax2.contourf(TT, PP, D12_masked * 1e9, levels=30, cmap='viridis')
cbar2 = fig.colorbar(cf2, ax=ax2, shrink=0.7)
cbar2.set_label(r'$D_{12}$ ($\times 10^{-9}$ m$^2$/s)', fontsize=9)
ax2.set_xlabel(r'$T$ (K)')
ax2.set_ylabel(r'$P$ (MPa)')
ax2.set_title(r'(b)  Diffusion $D_{12}(T,P)$')

# (c) mu contour
ax3 = fig.add_subplot(gs[1, 0])
cf3 = ax3.contourf(TT, PP, Mu_masked * 1000, levels=30, cmap='YlOrRd_r')
cbar3 = fig.colorbar(cf3, ax=ax3, shrink=0.7)
cbar3.set_label(r'$\mu_m$ (mPa$\cdot$s)', fontsize=9)
ax3.set_xlabel(r'$T$ (K)')
ax3.set_ylabel(r'$P$ (MPa)')
ax3.set_title(r'(c)  Viscosity $\mu_m(T,P)$')

# (d) Isothermal tau_ms vs P
ax4 = fig.add_subplot(gs[1, 1])
for i, T_sel in enumerate(iso_T):
    iT = np.argmin(np.abs(T_arr - T_sel))
    mask = ~np.isnan(Tau[iT, :])
    if np.sum(mask) < 2: continue
    ax4.plot(P_arr[mask], Tau[iT, mask], color=colors[i], lw=1.5, label=f'{T_sel:.0f} K')
ax4.set_xlabel(r'$P$ (MPa)')
ax4.set_ylabel(r'$\tau_{ms}$ (ms)')
ax4.set_title(r'(d)  $\tau_{ms}$ vs Pressure')
ax4.legend(fontsize=8, ncol=2)
ax4.grid(True, alpha=0.3)

# (e) Sensitivity
ax5 = fig.add_subplot(gs[1, 2])
ax5.loglog(D0_vals, tau_D0, 'b-', lw=1.5, label=r'$\tau_{ms} \propto 1/D_0$')
ax5.axvline(5e-8, color='red', ls='--', alpha=0.5, lw=0.8)
ax5.axhline(1.0, color='green', ls='--', alpha=0.5, lw=0.8)
ax5.axhline(100.0, color='orange', ls='--', alpha=0.5, lw=0.8)
ax5.set_xlabel(r'$D_0$ (m$^2$/s)')
ax5.set_ylabel(r'$\tau_{ms}$ (ms)')
ax5.set_title(r'(e)  Sensitivity to $D_0$')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3, which='both')

fig.suptitle(r'CO$_2$ + SW68:  Mesoscopic Transport — $\tau_{ms}$ Field  [PLACEHOLDER PARAMETERS]',
             fontsize=15, fontweight='bold', y=1.01)
fig.savefig('fig_taums_report.png', dpi=200)
plt.close()
print("[OK] fig_taums_report.png")

# ── Data table ────────────────────────────────────────────────────────────────
print("\n" + "=" * 90)
print("  Key tau_ms values (ms) — PLACEHOLDER PARAMETERS")
print("=" * 90)
header = f"{'T (K)':>8}"
for p in [3.0, 5.0, 8.0, 12.0, 15.0]:
    header += f"  P={p:.0f} MPa"
print(header)
print("-" * 90)
for T_sel in [283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15]:
    iT = np.argmin(np.abs(T_arr - T_sel))
    row = f"{T_sel:8.1f}"
    for P_sel in [3.0, 5.0, 8.0, 12.0, 15.0]:
        iP = np.argmin(np.abs(P_arr - P_sel))
        v = Tau[iT, iP]
        if np.isnan(v):
            row += "       full-misc"
        else:
            row += f"  {v:10.1f}"
    print(row)
print("-" * 90)
print(f"  tau_ms = L^2 / (D12 * Sh)  with L=10 mm, u=1.0 m/s")
print(f"  D0={5.0e-8:.1e} m2/s, gamma=0.6, V*_CO2={5.89e-4:.2e} m3/kg")
print(f"  To reach ~1 ms: increase D0 by ~1000x or decrease L by ~30x")
print("=" * 90)

print("\nAll figures saved:")
for f in ['fig_taums_contour.png', 'fig_d12_contour.png', 'fig_viscosity_contour.png',
           'fig_transport_panels.png', 'fig_sensitivity.png', 'fig_taums_report.png']:
    try:
        print(f"  {f:35s}  ({os.path.getsize(f)/1024:.0f} KB)")
    except:
        print(f"  {f:35s}  (MISSING)")
print("=" * 60)
