from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'figures'

d = pd.read_csv(ROOT / 'data' / 'processed' / 'PratxKapp_consistent_Dmax_scan.csv')
sm = pd.read_csv(ROOT / 'data' / 'processed' / 'PratxKapp_transfer_summary.csv')
conv = float(sm.loc[sm.design == 'CONV only', 'scaled_Jacobian_condition_number'].iloc[0])
bp = 4 / 0.42

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 26,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 19,
})

fig, ax = plt.subplots(figsize=(10.2, 5.8))
ax.plot(d.max_FLASH_dose_Gy, d.scaled_Jacobian_condition_number, marker='o', markersize=8.5, lw=3.2)
ax.axvline(bp, ls='--', lw=2.6, label=f'Depletion breakpoint\n({bp:.2f} Gy)')
ax.axhline(conv, ls=':', lw=2.6, label=f'CONV only\n({conv:.1f})')
ax.set_yscale('log')

ticks = [30, 100, 300, 600]
ax.yaxis.set_major_locator(FixedLocator(ticks))
ax.yaxis.set_major_formatter(FixedFormatter([str(x) for x in ticks]))

ax.set_xlabel('Maximum FLASH dose (Gy)')
ax.set_ylabel('Scaled Jacobian\ncondition number')
ax.grid(alpha=0.20, which='both')
ax.tick_params(axis='both', which='major', length=6, width=1.2)
ax.legend(frameon=True, framealpha=0.95, loc='upper right', bbox_to_anchor=(0.98, 0.82))

fig.tight_layout(pad=1.2)
fig.savefig(OUT / 'FigureS1.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT / 'FigureS1.pdf', bbox_inches='tight')
plt.close(fig)
