from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'figures';u=pd.read_csv(ROOT/'data'/'processed'/'prospective_uncertainty_summary_manuscript.csv')
fig,axes=plt.subplots(1,2,figsize=(13.5,4.8))
ax=axes[0];base=u[u.group=='survival_precision'].sort_values('survival_rel_SE');ax.plot(100*base.survival_rel_SE,base.half_width_pp,marker='o',lw=2.2);ax.set_xticks([1,3,5]);ax.set_ylim(0,7.5);ax.set_xlabel('Relative SE of clonogenic survival (%)');ax.set_ylabel('2.5–97.5% recovery half-width (pp)');ax.set_title('A. Effect of survival-measurement precision');ax.grid(alpha=.2)
ax=axes[1];d=u[u.group=='scenario'].copy();idx=np.arange(len(d));bars=ax.bar(idx,d.half_width_pp,width=.72)
# Highlight the key mode-differential-dose diagnostic.
for b,name in zip(bars,d.scenario):
    if name=='no_relative_dose': b.set_hatch('///'); b.set_linewidth(1.5)
labels=['Exact dose\nreference','Baseline dose\nuncertainty','Remove common\ndose error','Remove relative\ndose error','Remove oxygen\nerror','Remove survival\nnoise','Improved dose\ndelivery','Improved dose delivery\n+ O$_2$ precision','All errors\nzero']
ax.set_xticks(idx,labels,rotation=34,ha='right',fontsize=8.2);ax.set_ylabel('2.5–97.5% recovery half-width (pp)');ax.set_ylim(0,6.2);ax.set_title(r'B. Measurement-uncertainty scenarios ($\mathrm{SE}_{SF}=1\%$)')
for i,v in enumerate(d.half_width_pp):ax.text(i,v+.10,f'{v:.2f}',ha='center',fontsize=11,fontweight='bold')
fig.tight_layout();fig.savefig(OUT/'Figure4_final.png',dpi=300,bbox_inches='tight');fig.savefig(OUT/'Figure4_final.pdf',bbox_inches='tight');plt.close(fig)
