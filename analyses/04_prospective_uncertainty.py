import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from multiprocessing import Pool
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'processed'
OUT.mkdir(parents=True, exist_ok=True)
SEED=20260920
N_MC=int(os.environ.get('DC_RMLE_N_MC','1000'))
N_WORKERS=int(os.environ.get('DC_RMLE_WORKERS','6'))

# dc-RMLE reference / generator parameters
ALPHA1=0.006853; ALPHA2=0.166549; BETA=0.905239; KD=0.021876
GAMMA_C=0.130429; GAMMA_F=0.102507
K_O2=3.0; M_O2=3.0; P_NORM=142.6
O2_LEVELS=np.array([5.0,11.4,20.0])
DOSES=np.array([4.0,8.0,12.0,16.0])
G0_TRUE=0.4
LIPID_RESIDUAL_FRACTION=0.30

# Standard-uncertainty scenarios (k=1)
# u_rel is the SD of the FLASH/CONV relative dose ratio. Internally, independent
# mode offsets use u_mode = u_rel/sqrt(2), so the difference has SD u_rel.
SCENARIOS={
    'baseline': dict(u_common=0.020, u_rel=0.025, o2_scale=1.0),
    'improved_dose': dict(u_common=0.010, u_rel=0.010, o2_scale=1.0),
    'improved_all': dict(u_common=0.010, u_rel=0.010, o2_scale=0.5),
    'no_survival': dict(u_common=0.020, u_rel=0.025, o2_scale=1.0, sf_override=0.0),
    'no_common_dose': dict(u_common=0.0, u_rel=0.025, o2_scale=1.0),
    'no_relative_dose': dict(u_common=0.020, u_rel=0.0, o2_scale=1.0),
    'no_oxygen': dict(u_common=0.020, u_rel=0.025, o2_scale=0.0),
    'exact_dose_reference': dict(u_common=0.0, u_rel=0.0, o2_scale=1.0),
    'all_zero': dict(u_common=0.0, u_rel=0.0, o2_scale=0.0, sf_override=0.0),
}

# Pulse-structure constraint is a DESIGN restriction, not an uncertainty term:
# prospective use is restricted to an alanine-validated dose-per-pulse range.
ALANINE_DPP_VALIDATED_GY_PER_PULSE=(0.15,6.2)


def phi(p):
    p=np.maximum(np.asarray(p,float),0)
    return (p+K_O2/M_O2)/(p+K_O2)
Q_CONV=float(1/phi(P_NORM))


def o2_sd(p0):
    # Assumed baseline oxygen measurement performance; not a formal instrument spec.
    # 0.1 mmHg at 5 mmHg, increasing linearly toward 0.3 mmHg at 30 mmHg.
    if p0<=5: return 0.1
    return min(0.3,0.1+(p0-5)*(0.2/25.0))


def avg_phi(p0,D):
    frac=np.linspace(0,1,101)
    p=np.maximum(0,p0-G0_TRUE*D*frac)
    return float(np.trapezoid(phi(p),frac))


def m_from(D,lam):
    return (ALPHA1+ALPHA2*lam+BETA*lam*lam*(1-np.exp(-KD*D)))*D


def make_structure():
    rows=[]
    for state in ('vehicle','lipid_suppressed'):
        for p0 in O2_LEVELS:
            for mode in ('CONV','FLASH'):
                for D in DOSES:
                    rows.append((state,float(p0),mode,float(D)))
    return rows
STRUCT=make_structure()


def simulate_fit_one(args):
    scenario_name,sf_cv,true_delta,rep,seed=args
    cfg=SCENARIOS[scenario_name]
    sf_use=cfg.get('sf_override',sf_cv)
    rng=np.random.default_rng(seed)

    # Common random numbers: all latent draws are made for every scenario,
    # even when their scale is zero.
    z_common=rng.normal()
    z_mode={'CONV':rng.normal(),'FLASH':rng.normal()}
    z_cal={float(p):rng.normal() for p in O2_LEVELS}
    z_read_conv={(float(p),float(D),st):rng.normal() for st in ('vehicle','lipid_suppressed') for p in O2_LEVELS for D in DOSES}
    z_read_flash={(float(p),float(D),st):rng.normal(size=101) for st in ('vehicle','lipid_suppressed') for p in O2_LEVELS for D in DOSES}
    z_sf=rng.normal(size=len(STRUCT))

    u_mode=cfg['u_rel']/np.sqrt(2.0)
    eps_common=cfg['u_common']*z_common
    eps_mode={m:u_mode*z_mode[m] for m in ('CONV','FLASH')}
    qf_true=(1+true_delta)*Q_CONV

    states=[]; p0s=[]; modes=[]; Dnom=[]; aph_obs=[]; mobs=[]
    frac=np.linspace(0,1,101)
    sigm=np.sqrt(np.log(1+sf_use**2)) if sf_use>0 else 0.0

    for i,(state,p0,mode,D) in enumerate(STRUCT):
        f=1.0 if state=='vehicle' else LIPID_RESIDUAL_FRACTION
        gam=(GAMMA_C if mode=='CONV' else GAMMA_F)*f
        q=Q_CONV if mode=='CONV' else qf_true

        # Delivered-dose uncertainty affects the dose actually received; the fit
        # uses nominal/measured dose, mimicking unresolved delivery/transfer error.
        Dtrue=D*(1+eps_common+eps_mode[mode])
        Dtrue=max(Dtrue,1e-6)

        # True and measured oxygen modifier along the actual delivered trajectory.
        sig=o2_sd(p0)*cfg['o2_scale']
        cal=sig*z_cal[p0]
        if mode=='CONV':
            aph_t=float(phi(p0))
            p_meas=max(0,p0+cal+0.35*sig*z_read_conv[(p0,D,state)])
            aph_m=float(phi(p_meas))
        else:
            p_true=np.maximum(0,p0-G0_TRUE*Dtrue*frac)
            aph_t=float(np.trapezoid(phi(p_true),frac))
            p_meas=np.maximum(0,p_true+cal+0.35*sig*z_read_flash[(p0,D,state)])
            aph_m=float(np.trapezoid(phi(p_meas),frac))

        lam_t=q*aph_t+gam
        mtrue=m_from(Dtrue,lam_t)
        mobs_i=mtrue - sigm*z_sf[i] if sf_use>0 else mtrue

        states.append(state); p0s.append(p0); modes.append(mode); Dnom.append(D); aph_obs.append(aph_m); mobs.append(mobs_i)

    state=np.array(states); mode=np.array(modes); Dfit=np.array(Dnom,float); aph=np.array(aph_obs,float); y=np.array(mobs,float)
    isC=(mode=='CONV'); isV=(state=='vehicle')

    x0=np.array([ALPHA1,ALPHA2,BETA,KD,GAMMA_C,GAMMA_F,GAMMA_C*LIPID_RESIDUAL_FRACTION,GAMMA_F*LIPID_RESIDUAL_FRACTION,qf_true])
    lo=np.array([0,0,0,1e-5,0,0,0,0,0.3]); hi=np.array([1,1,10,1,2,2,2,2,2])

    def fun(x):
        a1,a2,b,kd,gcv,gfv,gcs,gfs,qf=x
        gam=np.where(isV,np.where(isC,gcv,gfv),np.where(isC,gcs,gfs))
        q=np.where(isC,Q_CONV,qf)
        lam=q*aph+gam
        pred=(a1+a2*lam+b*lam*lam*(1-np.exp(-kd*Dfit)))*Dfit
        # Keep objective well-defined in zero-survival-noise diagnostic.
        scale=max(np.sqrt(np.log(1+sf_cv**2)),1e-6)
        return (pred-y)/scale

    res=least_squares(fun,x0,bounds=(lo,hi),max_nfev=1500,xtol=1e-9,ftol=1e-9,gtol=1e-9)
    qhat=float(res.x[-1]); dhat=100*(qhat/Q_CONV-1)
    return scenario_name,sf_cv,true_delta,rep,dhat,2*res.cost,bool(res.success),int(res.active_mask[-1])


def summarize(df):
    out=[]
    for keys,g in df.groupby(['scenario','survival_rel_SE','true_delta_fraction'],sort=False):
        scen,sf,td=keys
        v=g.delta_hat_percent.to_numpy()
        q025,med,q975=np.percentile(v,[2.5,50,97.5])
        out.append(dict(scenario=scen,survival_rel_SE=sf,true_delta_percent=100*td,
                        mean_hat_percent=float(v.mean()),median_hat_percent=float(med),
                        bias_of_mean_pp=float(v.mean()-100*td),SD_hat_pp=float(v.std(ddof=1)),
                        q2_5_percent=float(q025),q97_5_percent=float(q975),
                        interval_width_pp=float(q975-q025),half_width_pp=float((q975-q025)/2),
                        interval_crosses_zero=bool(q025<=0<=q975),
                        optimizer_success_fraction=float(g.optimizer_success.mean()),
                        qf_bound_fraction=float((g.qf_active_mask!=0).mean()),N_MC=len(g)))
    return pd.DataFrame(out)


def main():
    jobs=[]
    # Main scenarios at 1% survival SE, truth -2%.
    for scen in SCENARIOS:
        for rep in range(N_MC):
            jobs.append((scen,0.01,-0.02,rep,SEED+rep))
    # Baseline survival-precision sensitivity (3%,5%).
    for sf in (0.03,0.05):
        for rep in range(N_MC):
            jobs.append(('baseline',sf,-0.02,rep,SEED+rep))
    # Truth-conditionality check at baseline 1% SE.
    for td in (-0.05,-0.10):
        for rep in range(N_MC):
            jobs.append(('baseline',0.01,td,rep,SEED+rep))

    with Pool(processes=N_WORKERS) as pool:
        results=list(pool.imap_unordered(simulate_fit_one,jobs,chunksize=8))
    df=pd.DataFrame(results,columns=['scenario','survival_rel_SE','true_delta_fraction','replicate','delta_hat_percent','weighted_SSE','optimizer_success','qf_active_mask'])
    df=df.sort_values(['scenario','survival_rel_SE','true_delta_fraction','replicate'])
    sm=summarize(df)

    df.to_csv(OUT/'prospective_uncertainty_samples.csv',index=False)
    sm.to_csv(OUT/'prospective_uncertainty_summary_full.csv',index=False)
    pd.DataFrame([
        dict(scenario=k,u_common_k1=v['u_common'],u_relative_CONV_FLASH_k1=v['u_rel'],u_mode_each_k1=v['u_rel']/np.sqrt(2),oxygen_uncertainty_multiplier=v['o2_scale'],survival_override=v.get('sf_override',np.nan))
        for k,v in SCENARIOS.items()
    ]).to_csv(OUT/'prospective_uncertainty_scenarios.csv',index=False)
    print(sm.to_string(index=False))

if __name__=='__main__':
    main()
