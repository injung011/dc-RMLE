"""Pratx–Kapp inverse-problem transfer check.

This is an inverse-problem use of the simplified oxygen-depletion formulation,
not a re-evaluation of the original forward-model study.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data'/'processed'; OUT.mkdir(parents=True,exist_ok=True)
THETA=np.array([0.15,1.63,0.26,0.42],float) # alpha0,R,phi,L_ROD
P_CONV=[40.,20.,4.]; DOSES=np.array([2.,4.,6.,8.,10.,12.,15.,20.]); P0_FLASH=4.

def oer(p,R,phi):
    p=np.maximum(np.asarray(p,float),0); return 1+R*(1-np.exp(-phi*p))
def alpha(p,th): return th[0]*oer(p,th[1],th[2])
def m_conv(D,p,th): return float(alpha(p,th)*D)
def m_flash(D,p0,th,n=4001):
    x=np.linspace(0,D,n); p=np.maximum(p0-th[3]*x,0); return float(np.trapezoid(alpha(p,th),x))
def jac(funcs,th=THETA):
    J=np.zeros((len(funcs),4)); th=np.asarray(th,float)
    for j in range(4):
        h=1e-6*max(abs(th[j]),1.0); a=th.copy();b=th.copy();a[j]+=h;b[j]-=h
        J[:,j]=[(f(a)-f(b))/(2*h) for f in funcs]
    return J
def cond(funcs):
    J=jac(funcs); norms=np.linalg.norm(J,axis=0); active=norms>1e-12; Js=J[:,active]/norms[active]; s=np.linalg.svd(Js,compute_uv=False); return float(s[0]/s[-1]),s

def base_conv_funcs():
    return [lambda th,p=p,D=D:m_conv(D,p,th) for p in P_CONV for D in DOSES]
def manuscript_designs():
    base=base_conv_funcs(); rows=[]; sv=[]
    configs=[('CONV only',[]),('CONV + FLASH pre-breakpoint',[2,4,6,8]),('CONV + FLASH post-breakpoint',DOSES.tolist())]
    for name,ds in configs:
        funcs=base+[lambda th,D=float(D):m_flash(D,P0_FLASH,th) for D in ds]
        c,s=cond(funcs); rows.append((name,len(funcs),c)); sv += [(name,i+1,float(x)) for i,x in enumerate(s)]
    funcs=base+[lambda th,D=float(D):m_flash(D,P0_FLASH,th) for D in DOSES]+[lambda th,D=float(D):m_conv(D,0,th) for D in DOSES]
    c,s=cond(funcs); rows.append(('CONV + FLASH + anoxic reference',len(funcs),c)); sv += [('CONV + FLASH + anoxic reference',i+1,float(x)) for i,x in enumerate(s)]
    return pd.DataFrame(rows,columns=['design','n_observations','scaled_Jacobian_condition_number']),pd.DataFrame(sv,columns=['design','singular_value_rank','singular_value'])

def dmax_scan():
    # One consistent Supplementary grid: eight equally spaced FLASH doses from 2 Gy to Dmax.
    vals=[4.,6.,8.,9.,9.5,9.52,10.,12.,15.,20.,25.]; base=base_conv_funcs(); rows=[]
    for dmax in vals:
        ds=np.linspace(2,dmax,8); funcs=base+[lambda th,D=float(D):m_flash(D,P0_FLASH,th) for D in ds]
        c,_=cond(funcs); rows.append((dmax,8,c))
    return pd.DataFrame(rows,columns=['max_FLASH_dose_Gy','n_FLASH_doses','scaled_Jacobian_condition_number'])

def high_o2_compensation():
    # Numerical fixed-alpha0 refit under saturated CONV conditions (40/60/80 mmHg).
    D=DOSES; design=[(p,float(d)) for p in (40.,60.,80.) for d in D]
    y=np.array([m_conv(d,p,THETA) for p,d in design])
    rows=[]
    for a0 in [0.10,0.15,0.20,0.30,0.39]:
        def fun(z):
            th=np.array([a0,z[0],z[1],THETA[3]])
            return np.array([m_conv(d,p,th) for p,d in design])-y
        r=least_squares(fun,[max(1e-6,THETA[0]*(1+THETA[1])/a0-1),THETA[2]],bounds=([0,1e-5],[10,2]),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=10000)
        R,phi=r.x; rows.append((a0,float(r.fun@r.fun),R,phi,a0*(1+R)))
    return pd.DataFrame(rows,columns=['fixed_alpha0_Gy-1','SSE','fitted_R','fitted_phi_mmHg-1','alpha0_times_1plusR'])

def main():
    sm,sv=manuscript_designs(); sm.to_csv(OUT/'PratxKapp_transfer_summary.csv',index=False); sv.to_csv(OUT/'PratxKapp_transfer_singular_values.csv',index=False)
    scan=dmax_scan(); scan.to_csv(OUT/'PratxKapp_consistent_Dmax_scan.csv',index=False)
    comp=high_o2_compensation(); comp.to_csv(OUT/'PratxKapp_alpha0_R_compensation.csv',index=False)
    print('breakpoint',P0_FLASH/THETA[3]); print(sm.to_string(index=False)); print(scan.to_string(index=False)); print(comp.to_string(index=False))
if __name__=='__main__': main()
