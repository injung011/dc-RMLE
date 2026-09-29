from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'figures';OUT.mkdir(parents=True,exist_ok=True)
prof=pd.read_csv(ROOT/'data'/'processed'/'deletion_profiles.csv');cont=pd.read_csv(ROOT/'data'/'processed'/'hypFLASH_qF_GammaF_contour_41x41.csv')

# Stack the three panels vertically so that each panel can use the full manuscript width.
fig,axes=plt.subplots(3,1,figsize=(8.6,13.8))

ax=axes[0]
for omit,label in [('norm-CONV',r'omit norm-CONV: fix $\Gamma_C$'),('norm-FLASH',r'omit norm-FLASH: fix $\Gamma_F$')]:
 d=prof[prof.omitted==omit];ax.plot(d.fixed_value,d.SSE,lw=2.4,label=label)
ax.set_xlabel(r'Fixed parameter ($\Gamma_C$ or $\Gamma_F$)')
ax.set_ylabel('SSE for remaining three conditions')
ax.set_title('A. Flat profile after condition deletion')
ax.legend(frameon=True,framealpha=1,fontsize=10.5,loc='center right')
ax.grid(alpha=.2)

ax=axes[1]
pv=cont.pivot(index='GammaF',columns='qFLASH',values='SSE').sort_index().sort_index(axis=1)
Q,G=np.meshgrid(pv.columns.to_numpy(float),pv.index.to_numpy(float))
Z=pv.to_numpy(float);zmin=float(np.nanmin(Z))
deltas=[.350,.120,.040,.010,.002]
styles=['solid',(0,(8,4)),(0,(6,3,1.5,3)),(0,(4,2)),(0,(2,1.5))]
for delta,ls in zip(deltas,styles):
 ax.contour(Q,G,Z,levels=[zmin+delta],linestyles=[ls],linewidths=2.2)
i=np.nanargmin(Z);ig,iq=np.unravel_index(i,Z.shape)
qb=pv.columns.to_numpy(float)[iq];gb=pv.index.to_numpy(float)[ig]
ax.scatter([qb],[gb],marker='*',s=190,facecolor='red',edgecolor='white',lw=1.2,zorder=6)
h=[Line2D([0],[0],color='black',linestyle=ls,lw=2,label=rf'$\Delta$SSE = {d:.3f}') for d,ls in zip(deltas,styles)]
h.append(Line2D([0],[0],marker='*',ls='none',markerfacecolor='red',markeredgecolor='white',markersize=13,label='Grid minimum'))
ax.legend(handles=h,frameon=True,framealpha=1,fontsize=10.3,loc='upper right')
ax.set_xlabel(r'$q_{\mathrm{FLASH}}$')
ax.set_ylabel(r'$\Gamma_F$')
ax.set_title('B. Compensation valley after omitting hypoxia–FLASH')
ax.grid(alpha=.12)

# Exact ridge panel from analytic rescaling.
ax=axes[2]
raw=pd.read_csv(ROOT/'data'/'raw'/'DU145_reconstructed.csv')
K=3.;MO=3.;PN=142.6;PH=11.4;G0=.4
phi=lambda p:(np.asarray(p)+K/MO)/(np.asarray(p)+K)
QC=float(1/phi(PN))
COND=raw.condition.to_numpy(str);D=raw.dose_Gy.to_numpy(float);Y=-np.log(raw.SF.to_numpy(float))
ISF=np.char.find(COND.astype(str),'FLASH')>=0;ISN=np.char.find(COND.astype(str),'norm')>=0
def avgphi(p0,d,f):
 if not f:return float(phi(p0))
 x=np.linspace(0,d,801);return float(np.trapezoid(phi(np.maximum(0,p0-G0*x)),x)/d)
APH=np.array([avgphi(PN if n else PH,d,f) for n,d,f in zip(ISN,D,ISF)])
XREF=np.array([.00685378,.16654861,.90524569,.02187576,.13042963,.10250709,.99178789])
LO=np.array([0,0,0,1e-6,0,0,.3]);HI=np.array([1,1,10,1,2,2,2])
def pred(x,mask):
 a1,a2,b,kd,gc,gf,qf=x
 dm=D[mask];am=APH[mask];fm=ISF[mask];nm=ISN[mask]
 q=np.where(fm,qf,QC);gam=np.where(nm,np.where(fm,gf,gc),0.)
 L=q*am+gam
 return (a1+a2*L+b*L**2*(1-np.exp(-kd*dm)))*dm
mask=COND!='hyp-CONV';held=COND=='hyp-CONV';gc0=.1304
free=[0,1,2,3,5,6];x0=XREF.copy();x0[4]=gc0
r=least_squares(lambda z:pred(np.array([z[0],z[1],z[2],z[3],gc0,z[4],z[5]]),mask)-Y[mask],
                x0[free],bounds=(LO[free],HI[free]),x_scale='jac',max_nfev=10000,
                xtol=1e-12,ftol=1e-12,gtol=1e-12)
xb=x0.copy();xb[free]=r.x
gcs=np.linspace(0,.99,500);sses=[];rm=[]
for gc in gcs:
 c=(1+gc)/(1+gc0)
 x=xb.copy();x[1]/=c;x[2]/=c**2;x[4]=gc;x[5]*=c;x[6]*=c
 sses.append(np.sum((pred(x,mask)-Y[mask])**2))
 rm.append(np.sqrt(np.mean((pred(x,held)-Y[held])**2)))
sses=np.array(sses);rm=np.array(rm)
imin=np.argmin(rm);gcmin=gcs[imin]
ax.plot(gcs,sses,lw=2.4,label='Three-condition fit SSE')
ax.axvline(gcmin,ls=':',lw=1.8,zorder=1)
ax.set_xlabel(r'Fixed $\Gamma_C$')
ax.set_ylabel('SSE for remaining three conditions')
ax.set_ylim(.2935,.2949)
ax.grid(alpha=.2)
ax2=ax.twinx()
ax2.plot(gcs,rm,ls='--',lw=2.4,label='Held-out hypoxia–CONV RMSE')
ax2.set_ylabel(r'Held-out RMSE in $m=-\ln(SF)$ space')
ax.set_title('C. Exact ridge and held-out prediction divergence')
ax.text(.03,.08,r'SSE range $<10^{-10}$',transform=ax.transAxes,fontsize=10.5)
h1,l1=ax.get_legend_handles_labels();h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,frameon=True,framealpha=1,fontsize=10.2,loc='upper left')
pd.DataFrame({'heldout_RMSE_min_GammaC':[gcmin],'heldout_RMSE_min':[rm[imin]],'SSE_range':[sses.max()-sses.min()]}).to_csv(ROOT/'data'/'processed'/'Figure2C_RMSE_minimum.csv',index=False)

fig.tight_layout(h_pad=2.4)
fig.savefig(OUT/'Figure2_final.png',dpi=300,bbox_inches='tight')
fig.savefig(OUT/'Figure2_final.pdf',bbox_inches='tight')
plt.close(fig)
