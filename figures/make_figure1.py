from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

du = pd.read_csv(ROOT / 'data' / 'raw' / 'DU145_reconstructed.csv')

K = 3.0
MO = 3.0
PN = 142.6
PH = 11.4
G0 = 0.4
phi = lambda p: (np.asarray(p) + K / MO) / (np.asarray(p) + K)
QC = float(1 / phi(PN))

a1, a2, b, kd, gc, gf, qf = [
    0.006853401311332844,
    0.16654888072896668,
    0.9052417566780149,
    0.021875890368583895,
    0.13042933779848265,
    0.10250685764143685,
    0.9917878827438733,
]


def av(p0, D):
    out = []
    for d in np.asarray(D, float):
        x = np.linspace(0, d, 1001)
        out.append(np.trapezoid(phi(np.maximum(0, p0 - G0 * x)), x) / d)
    return np.asarray(out)


def sf(c, D):
    D = np.asarray(D, float)
    if c == 'norm-CONV':
        L = np.full_like(D, 1 + gc)
    elif c == 'norm-FLASH':
        L = qf * av(PN, D) + gf
    elif c == 'hyp-CONV':
        L = np.full_like(D, QC * phi(PH))
    else:
        L = qf * av(PH, D)
    return np.exp(-(a1 + a2 * L + b * L**2 * (1 - np.exp(-kd * D))) * D)


lq = {}
rows = []
for c, d in du.groupby('condition', sort=False):
    D = d.dose_Gy.to_numpy(float)
    y = -np.log(d.SF.to_numpy(float))
    X = np.c_[D, D**2]
    coef, _ = nnls(X, y)
    lq[c] = coef
    rows.append((c, *coef, float(np.sum((X @ coef - y) ** 2))))

pd.DataFrame(
    rows,
    columns=['condition', 'alpha_LQ_Gy-1', 'beta_LQ_Gy-2', 'SSE_mspace']
).to_csv(ROOT / 'data' / 'processed' / 'LQ_benchmark.csv', index=False)

order = [
    ('hyp-CONV', '(A)'),
    ('hyp-FLASH', '(B)'),
    ('norm-CONV', '(C)'),
    ('norm-FLASH', '(D)'),
]

xlims = {
    'hyp-CONV': (0, 27),
    'hyp-FLASH': (0, 27),
    'norm-CONV': (0, 16.2),
    'norm-FLASH': (0, 16.2),
}

COLOR_DATA = '#1f77b4'
COLOR_DCRMLE = '#d62728'
COLOR_LQ = '#2ca02c'

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 26,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    # legend font size adjustment
    'legend.fontsize': 17,
})

fig, axs = plt.subplots(2, 2, figsize=(12.3, 10.2), sharey=True)

for ax, (c, panel_label) in zip(axs.ravel(), order):
    d = du[du.condition == c]
    xx = np.linspace(0.001, xlims[c][1], 450)
    err = d.graphical_error_SF.to_numpy(float)
    ok = np.isfinite(err)

    ax.errorbar(
        d.loc[ok, 'dose_Gy'], d.loc[ok, 'SF'], yerr=err[ok], fmt='o',
        ms=8.0, capsize=5.0, lw=1.8, color=COLOR_DATA, ecolor=COLOR_DATA,
        markerfacecolor=COLOR_DATA, markeredgecolor=COLOR_DATA,
        label='Recon. data', zorder=4
    )
    if np.any(~ok):
        ax.plot(d.loc[~ok, 'dose_Gy'], d.loc[~ok, 'SF'], 'o', ms=8.0, color=COLOR_DATA, zorder=4)

    ax.plot(xx, sf(c, xx), lw=3.2, color=COLOR_DCRMLE, label='dc-RMLE fit')
    al, bl = lq[c]
    ax.plot(xx, np.exp(-(al * xx + bl * xx**2)), lw=3.0, ls='--', color=COLOR_LQ, label='LQ fit')

    ax.set_yscale('log')
    ax.set_xlim(*xlims[c])
    ax.set_ylim(1e-5, 1.0)
    ax.set_xlabel('Dose (Gy)')
    ax.set_ylabel('Surviving fraction')
    ax.grid(alpha=0.18, which='both')
    ax.tick_params(axis='both', which='major', length=6, width=1.2)
    ax.tick_params(axis='y', labelleft=True)

    ax.text(
        0.5, 0.965, panel_label,
        transform=ax.transAxes,
        ha='center', va='top',
        fontsize=22, fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='none', pad=1.5)
    )

    if c in ('hyp-CONV', 'hyp-FLASH'):
        ax.set_xticks(np.arange(0, 28, 5))
    else:
        ax.set_xticks(np.arange(0, 17, 5))

    ax.legend(loc='lower left', frameon=True, framealpha=0.95)

fig.tight_layout(pad=1.4)
fig.savefig(OUT / 'Figure1.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT / 'Figure1.pdf', bbox_inches='tight')
plt.close(fig)
