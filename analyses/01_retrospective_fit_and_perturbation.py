"""Reference dc-RMLE fit and reported-only reconstruction perturbation.

The percentile ranges are reconstruction-sensitivity summaries under the stated
perturbation model, not formal confidence intervals.
"""
import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from multiprocessing import Pool

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "DU145_reconstructed.csv"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260921

K=3.0; MO=3.0; PN=142.6; PH=11.4; G0=0.4
phi=lambda p:(np.asarray(p)+K/MO)/(np.asarray(p)+K)
QC=float(1/phi(PN))
XREF=np.array([.00685378,.16654861,.90524569,.02187576,.13042963,.10250709,.99178789])
LO_BASE=np.array([0,0,0,1e-6,0,0,.3],float)


def prepare(beta_upper=10.0):
    df=pd.read_csv(RAW)
    cond=df.condition.to_numpy(str); D=df.dose_Gy.to_numpy(float); sf=df.SF.to_numpy(float); y=-np.log(sf)
    isf=np.char.find(cond.astype(str),'FLASH')>=0; isn=np.char.find(cond.astype(str),'norm')>=0
    def avgphi(p0,d,flash):
        if not flash:return float(phi(p0))
        x=np.linspace(0,d,501); p=np.maximum(0,p0-G0*x)
        return float(np.trapezoid(phi(p),x)/d)
    aph=np.array([avgphi(PN if n else PH,d,f) for n,d,f in zip(isn,D,isf)])
    hi=np.array([1,1,beta_upper,1,2,2,2],float)
    cv=np.array([0 if pd.isna(e) else e/s for e,s in zip(df.graphical_error_SF,sf)])
    sig=np.sqrt(np.log1p(cv**2))
    return df,cond,D,sf,y,isf,isn,aph,hi,sig


def fit_factory(beta_upper=10.0):
    df,cond,D,sf,y,isf,isn,aph,hi,sig=prepare(beta_upper)
    def model(x):
        a1,a2,b,kd,gc,gf,qf=x
        q=np.where(isf,qf,QC); gam=np.where(isn,np.where(isf,gf,gc),0.)
        lam=q*aph+gam
        return (a1+a2*lam+b*lam**2*(1-np.exp(-kd*D)))*D
    def fit_y(yy):
        starts=[np.minimum(XREF,hi-1e-10)]
        best=None
        for st in starts:
            r=least_squares(lambda x:model(x)-yy,st,bounds=(LO_BASE,hi),max_nfev=5000,
                            xtol=1e-11,ftol=1e-11,gtol=1e-11,x_scale='jac')
            if best is None or r.fun@r.fun < best.fun@best.fun: best=r
        return best
    return y,sig,hi,fit_y


def derive(x):
    R=(1+x[4])*float(phi(PN)/phi(PH))
    delta=100*(x[6]/QC-1)
    return R,delta


def worker(args):
    rep,beta_upper=args
    y,sig,hi,fit_y=fit_factory(beta_upper)
    rng=np.random.default_rng(SEED+rep)
    yy=y-sig*rng.normal(size=len(y))
    r=fit_y(yy); R,d=derive(r.x)
    return [rep,*r.x,R,d,float(r.fun@r.fun),bool(r.success)]


def run(n=1000,beta_upper=10.0,workers=6):
    y,sig,hi,fit_y=fit_factory(beta_upper)
    rc=fit_y(y); R0,d0=derive(rc.x)
    jobs=[(i,beta_upper) for i in range(n)]
    if workers==1: rows=[worker(x) for x in jobs]
    else:
        with Pool(workers) as p: rows=list(p.imap_unordered(worker,jobs,chunksize=5))
    cols=['rep','alpha1','alpha2','beta','kD','GammaC','GammaF','qFLASH','R_Lambda','Delta_percent','SSE','success']
    mc=pd.DataFrame(rows,columns=cols).sort_values('rep')
    suffix='' if beta_upper==10 else f'_beta{beta_upper:g}'
    mc.to_csv(OUT/f'retrospective_reported_only_samples{suffix}.csv',index=False)
    refs=dict(zip(cols[1:8],rc.x)); refs['R_Lambda']=R0; refs['Delta_percent']=d0
    summ=[]
    for col in cols[1:10]:
        v=mc[col].to_numpy(); q=np.percentile(v,[2.5,50,97.5])
        summ.append((col,refs.get(col,np.nan),q[0],q[1],q[2],float(v.mean()),float(v.std(ddof=1))))
    sm=pd.DataFrame(summ,columns=['quantity','reference','q2_5','median','q97_5','mean','sd'])
    sm.to_csv(OUT/f'retrospective_reported_only_summary{suffix}.csv',index=False)
    print('central SSE',float(rc.fun@rc.fun),'R_Lambda',R0,'Delta_percent',d0)
    print(sm.to_string(index=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,default=1000); ap.add_argument('--beta-upper',type=float,default=10.0); ap.add_argument('--workers',type=int,default=6)
    a=ap.parse_args(); run(a.n,a.beta_upper,a.workers)
