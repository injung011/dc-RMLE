from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'processed'
OUT.mkdir(parents=True, exist_ok=True)

RAW = [
('norm-CONV',3.05,.505,.2),('norm-CONV',6.70,.09,.07),('norm-CONV',9.82,.0132,.007),('norm-CONV',12.78,.0022,.00205),('norm-CONV',14.98,.00048,.000313),
('norm-FLASH',.62,.91,None),('norm-FLASH',3.50,.43,None),('norm-FLASH',6.64,.11,.04),('norm-FLASH',9.39,.028,.017),('norm-FLASH',13.07,.0029,.0021),('norm-FLASH',15.43,.000505,.00048),
('hyp-CONV',6.67,.167,.14),('hyp-CONV',9.62,.058,.03),('hyp-CONV',12.89,.0179,.015),('hyp-CONV',15.23,.006,.004),('hyp-CONV',17.93,.000896,.0004),('hyp-CONV',21.20,.000193,.000128),
('hyp-FLASH',.60,1.0,.2),('hyp-FLASH',3.55,.573,.28),('hyp-FLASH',6.17,.224,.155),('hyp-FLASH',9.55,.098,.08),('hyp-FLASH',12.94,.02,.0127),('hyp-FLASH',15.34,.00838,.006),('hyp-FLASH',18.01,.00342,.0019),('hyp-FLASH',21.34,.00051,.0003),('hyp-FLASH',25,.000157,7.5e-5)]
DF = pd.DataFrame(RAW, columns=['condition','dose','SF','graph_err'])
K=3.0; MO=3.0; PN=142.6; PH=11.4; G0=0.4
phi=lambda p:(np.asarray(p)+K/MO)/(np.asarray(p)+K)
QC=float(1/phi(PN))

def avgphi(p0,D,flash):
    if not flash: return float(phi(p0))
    x=np.linspace(0,D,801); p=np.maximum(0,p0-G0*x)
    return float(np.trapezoid(phi(p),x)/D)

COND=DF.condition.to_numpy(); D=DF.dose.to_numpy(float); SF=DF.SF.to_numpy(float); Y=-np.log(SF)
ISF=np.char.find(COND.astype(str),'FLASH')>=0; ISN=np.char.find(COND.astype(str),'norm')>=0
APH=np.array([avgphi(PN if n else PH,d,f) for n,d,f in zip(ISN,D,ISF)])
XREF=np.array([.00685378,.16654861,.90524569,.02187576,.13042963,.10250709,.99178789])
LO=np.array([0,0,0,1e-6,0,0,.3]); HI=np.array([1,1,10,1,2,2,2])

def predict(x, mask):
    a1,a2,b,kd,gc,gf,qf=x
    dm=D[mask]; am=APH[mask]; fm=ISF[mask]; nm=ISN[mask]
    q=np.where(fm,qf,QC); gam=np.where(nm,np.where(fm,gf,gc),0.)
    lam=q*am+gam
    return (a1+a2*lam+b*lam**2*(1-np.exp(-kd*dm)))*dm

def profile_fixed(omit, idx, grid):
    mask=COND!=omit; free=[i for i in range(7) if i!=idx]
    rows=[]
    for val in grid:
        x0=XREF.copy(); x0[idx]=val
        def fun(z):
            x=x0.copy(); x[free]=z
            return predict(x,mask)-Y[mask]
        r=least_squares(fun,x0[free],bounds=(LO[free],HI[free]),max_nfev=5000,xtol=1e-11,ftol=1e-11,gtol=1e-11,x_scale='jac')
        x=x0.copy(); x[free]=r.x
        rows.append((omit,val,float(r.fun@r.fun),*x))
    return rows

def contour(nq=41, ng=41):
    # Deterministic profile-SSE surface. No random seed is involved.
    mask=COND!='hyp-FLASH'; free=[0,1,2,3,4]
    qs=np.linspace(0.60,1.32,nq); gs=np.linspace(0.0,0.55,ng)
    rows=[]
    for q in qs:
        warm=XREF[free].copy()
        for g in gs:
            def fun(z):
                x=np.array([z[0],z[1],z[2],z[3],z[4],g,q])
                return predict(x,mask)-Y[mask]
            r=least_squares(fun,warm,bounds=(LO[free],HI[free]),max_nfev=700,xtol=1e-7,ftol=1e-7,gtol=1e-7,x_scale='jac')
            warm=r.x
            rows.append((q,g,float(r.fun@r.fun)))
    return pd.DataFrame(rows,columns=['qFLASH','GammaF','SSE'])

def exact_ridge_table():
    # Fit the reduced design once at the manuscript reference Gamma_C, then
    # generate the complete exact ridge using the analytic reduced-design rescaling.
    mask=COND!='hyp-CONV'; gc0=0.1304; free=[0,1,2,3,5,6]
    x0=XREF.copy(); x0[4]=gc0
    def fun(z):
        x=x0.copy(); x[free]=z
        return predict(x,mask)-Y[mask]
    r=least_squares(fun,x0[free],bounds=(LO[free],HI[free]),max_nfev=10000,
                    xtol=1e-12,ftol=1e-12,gtol=1e-12,x_scale='jac')
    xb=x0.copy(); xb[free]=r.x
    held=(COND=='hyp-CONV')
    rows=[]
    for gc in [0.0000,0.0500,0.1000,0.1304,0.1835,0.3000,0.6490,0.9900]:
        c=(1+gc)/(1+gc0)
        x=xb.copy()
        x[1]=xb[1]/c
        x[2]=xb[2]/c**2
        x[4]=gc
        x[5]=c*xb[5]
        x[6]=c*xb[6]
        sse=float(np.sum((predict(x,mask)-Y[mask])**2))
        R=(1+gc)*float(phi(PN)/phi(PH))
        delta=100*(x[6]/QC-1)
        rmse=float(np.sqrt(np.mean((predict(x,held)-Y[held])**2)))
        rows.append((gc,sse,R,delta,rmse))
    return pd.DataFrame(rows,columns=['fixed_GammaC','three_condition_SSE','R_Lambda',
        'DeltaPs_over_Ps_percent','heldout_hypCONV_RMSE_mspace'])

def main():
    prof=[]
    prof += profile_fixed('norm-CONV',4,np.linspace(0,1.5,31))
    prof += profile_fixed('norm-FLASH',5,np.linspace(0,1.5,31))
    pd.DataFrame(prof,columns=['omitted','fixed_value','SSE','alpha1','alpha2','beta','kD','GammaC','GammaF','qFLASH']).to_csv(OUT/'deletion_profiles.csv',index=False)
    cont=contour(); cont.to_csv(OUT/'hypFLASH_qF_GammaF_contour_41x41.csv',index=False)
    exact_ridge_table().to_csv(OUT/'Table3_ridge_points.csv',index=False)
    print('contour points',len(cont),'SSE min',cont.SSE.min(),'max',cont.SSE.max())

if __name__=='__main__':
    main()
