from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'figures'

u = pd.read_csv(ROOT / 'data' / 'processed' / 'prospective_uncertainty_summary_manuscript.csv')

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 26,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 19,
})

fig, axes = plt.subplots(2, 1, figsize=(11.0, 12.8))

# ---------------- Panel A: scenarios ----------------
ax = axes[0]
d = u[u.group == 'scenario'].copy()
idx = np.arange(len(d))
bars = ax.bar(idx, d.half_width_pp, width=0.72)

for b, name in zip(bars, d.scenario):
    if name == 'no_relative_dose':
        b.set_hatch('///')
        b.set_linewidth(1.8)

labels = [
    'Exact dose ref.',
    'Baseline',
    'No common dose',
    'No rel. dose',
    'No O$_2$ error',
    'No surv. noise',
    'Improved dose',
    'Improved dose + O$_2$',
    'Zero-error',
]

ax.set_xticks(idx, labels, rotation=28, ha='right', rotation_mode='anchor')
ax.set_ylabel('2.5–97.5% recovery\nhalf-width (pp)')
ax.set_ylim(0, 6.4)
ax.grid(alpha=0.18, axis='y')
ax.tick_params(axis='x', pad=10)
ax.tick_params(axis='both', which='major', length=6, width=1.2)

for i, v in enumerate(d.half_width_pp):
    ax.text(i, v + 0.12, f'{v:.2f}', ha='center', fontsize=22, fontweight='bold')

ax.text(
    0.5, 0.965, '(A)', transform=ax.transAxes,
    ha='center', va='top', fontsize=22, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
)

# ---------------- Panel B: survival precision ----------------
ax = axes[1]
base = u[u.group == 'survival_precision'].sort_values('survival_rel_SE')
ax.plot(100 * base.survival_rel_SE, base.half_width_pp, marker='o', markersize=8.5, lw=3.2)
ax.set_xticks([1, 3, 5])
ax.set_ylim(0, 7.5)
ax.set_xlabel('Relative SE of clonogenic survival (%)')
ax.set_ylabel('2.5–97.5% recovery\nhalf-width (pp)')
ax.grid(alpha=0.20)
ax.tick_params(axis='both', which='major', length=6, width=1.2)
ax.text(
    0.5, 0.965, '(B)', transform=ax.transAxes,
    ha='center', va='top', fontsize=22, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
)

fig.tight_layout(pad=1.5)
fig.savefig(OUT / 'Figure4_final.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT / 'Figure4_final.pdf', bbox_inches='tight')
plt.close(fig)
