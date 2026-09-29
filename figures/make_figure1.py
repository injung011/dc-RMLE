from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import nnls
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'/'figures'; OUT.mkdir(parents=True,exist_ok=True)
du=pd.read_csv(ROOT/'data'/'raw'/'DU145_reconstructed.csv')
K=3.;MO=3.;PN=142.6;PH=11.4;G0=.4;phi=lambda p:(np.asarray(p)+K/MO)/(np.asarray(p)+K);QC=float(1/phi(PN))
a1,a2,b,kd,gc,gf,qf=[.006853401311332844,.16654888072896668,.9052417566780149,.021875890368583895,.13042933779848265,.10250685764143685,.9917878827438733]
def av(p0,D):
    out=[]
    for d in np.asarray(D,float):
        x=np.linspace(0,d,1001); out.append(np.trapezoid(phi(np.maximum(0,p0-G0*x)),x)/d)
    return np.asarray(out)
def sf(c,D):
    D=np.asarray(D,float)
    if c=='norm-CONV': L=np.full_like(D,1+gc)
    elif c=='norm-FLASH': L=qf*av(PN,D)+gf
    elif c=='hyp-CONV': L=np.full_like(D,QC*phi(PH))
    else: L=qf*av(PH,D)
    return np.exp(-(a1+a2*L+b*L**2*(1-np.exp(-kd*D)))*D)
lq={}; rows=[]
for c,d in du.groupby('condition',sort=False):
    D=d.dose_Gy.to_numpy(float);y=-np.log(d.SF.to_numpy(float));X=np.c_[D,D**2];coef,_=nnls(X,y);lq[c]=coef;rows.append((c,*coef,float(np.sum((X@coef-y)**2))))
pd.DataFrame(rows,columns=['condition','alpha_LQ_Gy-1','beta_LQ_Gy-2','SSE_mspace']).to_csv(ROOT/'data'/'processed'/'LQ_benchmark.csv',index=False)
order=[('hyp-CONV','A. Hypoxia–CONV'),('hyp-FLASH','B. Hypoxia–FLASH'),('norm-CONV','C. Normoxia–CONV'),('norm-FLASH','D. Normoxia–FLASH')]
xlims={'hyp-CONV':(0,26),'hyp-FLASH':(0,26),'norm-CONV':(0,16.2),'norm-FLASH':(0,16.2)}
fig,axs=plt.subplots(2,2,figsize=(10.4,8.4),sharey=True)
for ax,(c,title) in zip(axs.ravel(),order):
    d=du[du.condition==c]; xx=np.linspace(.001,xlims[c][1],450);err=d.graphical_error_SF.to_numpy(float);ok=np.isfinite(err)
    ax.errorbar(d.loc[ok,'dose_Gy'],d.loc[ok,'SF'],yerr=err[ok],fmt='o',capsize=3,label='Reconstructed data',zorder=4)
    if np.any(~ok): ax.plot(d.loc[~ok,'dose_Gy'],d.loc[~ok,'SF'],'o',zorder=4)
    ax.plot(xx,sf(c,xx),lw=2.2,label='dc-RMLE fit'); al,bl=lq[c]; ax.plot(xx,np.exp(-(al*xx+bl*xx**2)),lw=2,ls='--',label='LQ fit')
    ax.set_yscale('log');ax.set_xlim(*xlims[c]);ax.set_ylim(8e-5,1.5);ax.set_title(title);ax.set_xlabel('Dose (Gy)');ax.set_ylabel('Surviving fraction');ax.grid(alpha=.18,which='both');ax.legend(frameon=False,fontsize=9)
fig.tight_layout();fig.savefig(OUT/'Figure1_final.png',dpi=300,bbox_inches='tight');fig.savefig(OUT/'Figure1_final.pdf',bbox_inches='tight');plt.close(fig)
