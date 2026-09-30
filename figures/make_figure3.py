from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'figures'

d = pd.read_csv(ROOT / 'data' / 'processed' / 'oxygen_trajectory_sensitivity.csv')

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 26,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 17,
})

fig, ax = plt.subplots(figsize=(8.0, 5.8))
ax.plot(d.g0_mmHg_per_Gy, d.DeltaPs_over_Ps_percent, marker='o', markersize=8.5, lw=3.2)

ax.axvspan(0.16, 0.17, alpha=0.18, label='Electron FLASH\nmeasurement (0.16–0.17)')
ax.axvline(0.37, ls='--', lw=2.6, label='Proton FLASH\nmeasurement (~0.37)')
ax.axhline(0, lw=1.6, ls='--')

ax.set_xlabel('Prescribed FLASH oxygen-depletion rate\n'
              r'$g_0$ (mmHg Gy$^{-1}$)')
ax.set_ylabel(r'Refitted $\Delta P_s/P_s$ (%)')
ax.grid(alpha=0.20)
ax.tick_params(axis='both', which='major', length=6, width=1.2)
ax.legend(frameon=True, framealpha=0.95, loc='upper left', bbox_to_anchor=(0.02, 0.92))

fig.tight_layout(pad=1.2)
fig.savefig(OUT / 'Figure3_final.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT / 'Figure3_final.pdf', bbox_inches='tight')
plt.close(fig)
