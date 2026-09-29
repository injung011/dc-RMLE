from pathlib import Path
import pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'figures';d=pd.read_csv(ROOT/'data'/'processed'/'oxygen_trajectory_sensitivity.csv')
fig,ax=plt.subplots(figsize=(6.5,4.6));ax.plot(d.g0_mmHg_per_Gy,d.DeltaPs_over_Ps_percent,marker='o',lw=2.2)
for x,y in zip(d.g0_mmHg_per_Gy,d.DeltaPs_over_Ps_percent):ax.text(x,y+.35,f'{y:.2f}%',ha='center',fontsize=9)
# Literature-measured ranges shown as contextual references, not fitted intervals.
ax.axvspan(.16,.17,alpha=.18,label='Electron FLASH measurement (0.16–0.17)');ax.axvline(.37,ls='--',lw=1.5,label='Proton FLASH measurement (~0.37)')
ax.axhline(0,lw=.9,ls='--');ax.set_xlabel(r'Prescribed FLASH oxygen-depletion rate $g_0$ (mmHg Gy$^{-1}$)');ax.set_ylabel(r'Refitted $\Delta P_s/P_s$ (%)');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8.5,loc='lower right');fig.tight_layout();fig.savefig(OUT/'Figure3_final.png',dpi=300,bbox_inches='tight');fig.savefig(OUT/'Figure3_final.pdf',bbox_inches='tight');plt.close(fig)
