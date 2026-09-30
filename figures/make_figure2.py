from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

prof = pd.read_csv(ROOT / 'data' / 'processed' / 'deletion_profiles.csv')
cont = pd.read_csv(ROOT / 'data' / 'processed' / 'hypFLASH_qF_GammaF_contour_41x41.csv')

plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 24,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 16,
})

fig, axes = plt.subplots(3, 1, figsize=(8.6, 12.6))

# ---------------- Panel A ----------------
ax = axes[0]
for omit, label in [
    ('norm-CONV', r'omit norm-CONV: fix $\Gamma_C$'),
    ('norm-FLASH', r'omit norm-FLASH: fix $\Gamma_F$')
]:
    d = prof[prof.omitted == omit]
    ax.plot(d.fixed_value, d.SSE, lw=2.8, label=label)

ax.set_xlabel(r'Fixed parameter ($\Gamma_C$ or $\Gamma_F$)')
ax.set_ylabel('SSE for remaining\nthree conditions')
ax.set_ylim(0.30, 0.50)
ax.set_yticks(np.arange(0.30, 0.51, 0.05))
ax.grid(alpha=0.20)
ax.tick_params(axis='both', which='major', length=6, width=1.2)
ax.legend(frameon=True, framealpha=0.95, loc='lower right')
ax.text(
    0.5, 0.965, '(A)', transform=ax.transAxes,
    ha='center', va='top', fontsize=22, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
)

# ---------------- Panel B ----------------
ax = axes[1]
pv = cont.pivot(index='GammaF', columns='qFLASH', values='SSE').sort_index().sort_index(axis=1)
Q, G = np.meshgrid(pv.columns.to_numpy(float), pv.index.to_numpy(float))
Z = pv.to_numpy(float)
zmin = float(np.nanmin(Z))
deltas = [0.350, 0.120, 0.040, 0.010]
styles = ['solid', (0, (8, 4)), (0, (6, 3, 1.5, 3)), (0, (4, 2))]

for delta, ls in zip(deltas, styles):
    ax.contour(Q, G, Z, levels=[zmin + delta], linestyles=[ls], linewidths=2.4)

i = np.nanargmin(Z)
ig, iq = np.unravel_index(i, Z.shape)
qb = pv.columns.to_numpy(float)[iq]
gb = pv.index.to_numpy(float)[ig]
ax.scatter([qb], [gb], marker='*', s=440, facecolor='red', edgecolor='white', lw=1.2, zorder=6)

handles = [
    Line2D([0], [0], color='black', linestyle=ls, lw=2.4, label=rf'$\Delta$SSE = {d:.3f}')
    for d, ls in zip(deltas, styles)
]
handles.append(
    Line2D([0], [0], marker='*', ls='none', markerfacecolor='red', markeredgecolor='white',
           markersize=14, label='Grid minimum')
)

legend_main = ax.legend(handles=handles[:-1], frameon=True, framealpha=0.95, loc='upper right', fontsize=14)
ax.add_artist(legend_main)
ax.legend(handles=[handles[-1]], frameon=True, framealpha=0.95, loc='lower left', fontsize=14)
ax.set_xlabel(r'$q_{\mathrm{FLASH}}$')
ax.set_ylabel(r'$\Gamma_F$')
ax.set_ylim(0, 0.6)
ax.set_yticks(np.arange(0.0, 0.61, 0.1))
ax.grid(alpha=0.12)
ax.tick_params(axis='both', which='major', length=6, width=1.2)
ax.text(
    0.5, 0.965, '(B)', transform=ax.transAxes,
    ha='center', va='top', fontsize=22, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
)

# ---------------- Panel C ----------------
ax = axes[2]
raw = pd.read_csv(ROOT / 'data' / 'raw' / 'DU145_reconstructed.csv')
K = 3.0
MO = 3.0
PN = 142.6
PH = 11.4
G0 = 0.4
phi = lambda p: (np.asarray(p) + K / MO) / (np.asarray(p) + K)
QC = float(1 / phi(PN))
COND = raw.condition.to_numpy(str)
D = raw.dose_Gy.to_numpy(float)
Y = -np.log(raw.SF.to_numpy(float))
ISF = np.char.find(COND.astype(str), 'FLASH') >= 0
ISN = np.char.find(COND.astype(str), 'norm') >= 0


def avgphi(p0, d, f):
    if not f:
        return float(phi(p0))
    x = np.linspace(0, d, 801)
    return float(np.trapezoid(phi(np.maximum(0, p0 - G0 * x)), x) / d)

APH = np.array([avgphi(PN if n else PH, d, f) for n, d, f in zip(ISN, D, ISF)])
XREF = np.array([0.00685378, 0.16654861, 0.90524569, 0.02187576, 0.13042963, 0.10250709, 0.99178789])
LO = np.array([0, 0, 0, 1e-6, 0, 0, 0.3])
HI = np.array([1, 1, 10, 1, 2, 2, 2])


def pred(x, mask):
    a1, a2, b, kd, gc, gf, qf = x
    dm = D[mask]
    am = APH[mask]
    fm = ISF[mask]
    nm = ISN[mask]
    q = np.where(fm, qf, QC)
    gam = np.where(nm, np.where(fm, gf, gc), 0.0)
    L = q * am + gam
    return (a1 + a2 * L + b * L**2 * (1 - np.exp(-kd * dm))) * dm

mask = COND != 'hyp-CONV'
held = COND == 'hyp-CONV'
gc0 = 0.1304
free = [0, 1, 2, 3, 5, 6]
x0 = XREF.copy()
x0[4] = gc0
r = least_squares(
    lambda z: pred(np.array([z[0], z[1], z[2], z[3], gc0, z[4], z[5]]), mask) - Y[mask],
    x0[free], bounds=(LO[free], HI[free]), x_scale='jac',
    max_nfev=10000, xtol=1e-12, ftol=1e-12, gtol=1e-12
)
xb = x0.copy()
xb[free] = r.x

gcs = np.linspace(0, 0.99, 500)
sses = []
rm = []
for gc in gcs:
    c = (1 + gc) / (1 + gc0)
    x = xb.copy()
    x[1] /= c
    x[2] /= c**2
    x[4] = gc
    x[5] *= c
    x[6] *= c
    sses.append(np.sum((pred(x, mask) - Y[mask]) ** 2))
    rm.append(np.sqrt(np.mean((pred(x, held) - Y[held]) ** 2)))

sses = np.array(sses)
rm = np.array(rm)
imin = np.argmin(rm)
gcmin = gcs[imin]

ax.plot(gcs, sses, lw=2.8, label='Three-condition fit SSE')
ax.axvline(gcmin, ls=':', lw=2.0, zorder=1)
ax.set_xlabel(r'Fixed $\Gamma_C$')
ax.set_ylabel('SSE for remaining\nthree conditions')
ax.set_ylim(0.2935, 0.2950)
ax.set_yticks(np.arange(0.2935, 0.2950 + 1e-9, 0.0005))
ax.grid(alpha=0.20)
ax.tick_params(axis='both', which='major', length=6, width=1.2)

ax2 = ax.twinx()
ax2.plot(gcs, rm, ls='--', lw=2.8, label='Held-out hypoxia–CONV RMSE')
ax2.set_ylim(bottom=0)
ax2.set_ylabel('Held-out RMSE\nin ' + r'$m=-\ln(SF)$ space', fontsize=20)
ax2.tick_params(axis='y', labelsize=18)

ax.text(
    0.5, 0.965, '(C)', transform=ax.transAxes,
    ha='center', va='top', fontsize=22, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=True, framealpha=0.95, fontsize=14, loc='lower right')

pd.DataFrame({
    'heldout_RMSE_min_GammaC': [gcmin],
    'heldout_RMSE_min': [rm[imin]],
    'SSE_range': [sses.max() - sses.min()]
}).to_csv(ROOT / 'data' / 'processed' / 'Figure2C_RMSE_minimum.csv', index=False)

fig.tight_layout(pad=1.4)
fig.savefig(OUT / 'Figure2.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT / 'Figure2.pdf', bbox_inches='tight')
plt.close(fig)
