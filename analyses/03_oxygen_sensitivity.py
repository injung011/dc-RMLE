from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data'/'processed'; OUT.mkdir(parents=True,exist_ok=True)
raw=pd.read_csv(ROOT/'data'/'raw'/'DU145_reconstructed.csv')
K=3.; MO=3.; PN=142.6; PH=11.4
phi=lambda p:(np.asarray(p)+K/MO)/(np.asarray(p)+K); QC=float(1/phi(PN))
X0=np.array([.00685378,.16654861,.90524569,.02187576,.13042963,.10250709,.99178789])
LO=np.array([0,0,0,1e-6,0,0,.3]); HI=np.array([1,1,10,1,2,2,2])
cond=raw.condition.to_numpy(); D=raw.dose_Gy.to_numpy(float); Y=-np.log(raw.SF.to_numpy(float)); isf=np.char.find(cond.astype(str),'FLASH')>=0; isn=np.char.find(cond.astype(str),'norm')>=0

def fit(g0):
    def avgphi(p0,d,flash):
        if not flash:return float(phi(p0))
        x=np.linspace(0,d,801); p=np.maximum(0,p0-g0*x); return float(np.trapezoid(phi(p),x)/d)
    aph=np.array([avgphi(PN if n else PH,d,f) for n,d,f in zip(isn,D,isf)])
    def pred(x):
        a1,a2,b,kd,gc,gf,qf=x; q=np.where(isf,qf,QC); gam=np.where(isn,np.where(isf,gf,gc),0.); lam=q*aph+gam
        return (a1+a2*lam+b*lam*lam*(1-np.exp(-kd*D)))*D
    r=least_squares(lambda x:pred(x)-Y,X0,bounds=(LO,HI),x_scale='jac',max_nfev=10000,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    delta=100*(r.x[6]/QC-1)
    return [g0,delta,float(r.fun@r.fun),*r.x]
rows=[fit(g) for g in [0.10,0.20,0.40]]
pd.DataFrame(rows,columns=['g0_mmHg_per_Gy','DeltaPs_over_Ps_percent','SSE','alpha1','alpha2','beta','kD','GammaC','GammaF','qFLASH']).to_csv(OUT/'oxygen_trajectory_sensitivity.csv',index=False)
