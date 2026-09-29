from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data'/'processed'
rows=[]
for bound,suffix in [(5.0,'_beta5'),(10.0,'')]:
    p=OUT/f'retrospective_reported_only_samples{suffix}.csv'
    if not p.exists(): raise FileNotFoundError(p)
    d=pd.read_csv(p); rq=np.percentile(d.R_Lambda,[2.5,50,97.5]); dq=np.percentile(d.Delta_percent,[2.5,50,97.5])
    rows.append(dict(beta_upper_bound_Gy_1=bound,beta_bound_occupancy_percent=100*np.mean(d.beta>=bound-1e-6),R_Lambda_q2_5=rq[0],R_Lambda_median=rq[1],R_Lambda_q97_5=rq[2],Delta_percent_q2_5=dq[0],Delta_percent_median=dq[1],Delta_percent_q97_5=dq[2]))
pd.DataFrame(rows).to_csv(OUT/'beta_bound_sensitivity.csv',index=False)
