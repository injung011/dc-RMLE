from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'figures'
u=pd.read_csv(ROOT/'data'/'processed'/'prospective_uncertainty_summary_manuscript.csv')

# Stack panels vertically for readability. Scenario comparison is panel A so the
# figure follows the order used in the main Results narrative.
fig,axes=plt.subplots(2,1,figsize=(9.0,10.8))

ax=axes[0]
d=u[u.group=='scenario'].copy()
idx=np.arange(len(d))
bars=ax.bar(idx,d.half_width_pp,width=.72)
for b,name in zip(bars,d.scenario):
    if name=='no_relative_dose':
        b.set_hatch('///'); b.set_linewidth(1.5)
labels=[
    'Exact dose\nreference',
    'Baseline dose\nuncertainty',
    'Remove common\ndose error',
    'Remove relative\ndose error',
    'Remove oxygen\nerror',
    'Remove survival\nnoise',
    'Improved dose\ndelivery',
    'Improved dose delivery\n+ O$_2$ precision',
    'All errors\nzero'
]
ax.set_xticks(idx,labels,rotation=28,ha='right',fontsize=9.2)
ax.set_ylabel('2.5–97.5% recovery half-width (pp)')
ax.set_ylim(0,6.2)
ax.set_title(r'A. Measurement-uncertainty scenarios ($\mathrm{SE}_{SF}=1\%$)')
for i,v in enumerate(d.half_width_pp):
    ax.text(i,v+.10,f'{v:.2f}',ha='center',fontsize=10.5,fontweight='bold')
ax.grid(axis='y',alpha=.2)

ax=axes[1]
base=u[u.group=='survival_precision'].sort_values('survival_rel_SE')
ax.plot(100*base.survival_rel_SE,base.half_width_pp,marker='o',lw=2.4)
for x,y in zip(100*base.survival_rel_SE,base.half_width_pp):
    ax.text(x,y+.16,f'{y:.2f}',ha='center',fontsize=10.5)
ax.set_xticks([1,3,5])
ax.set_ylim(0,7.7)
ax.set_xlabel('Relative SE of clonogenic survival (%)')
ax.set_ylabel('2.5–97.5% recovery half-width (pp)')
ax.set_title('B. Effect of survival-measurement precision')
ax.grid(alpha=.2)

fig.tight_layout(h_pad=2.6)
fig.savefig(OUT/'Figure4_final.png',dpi=300,bbox_inches='tight')
fig.savefig(OUT/'Figure4_final.pdf',bbox_inches='tight')
plt.close(fig)
