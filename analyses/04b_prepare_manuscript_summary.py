from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data'/'processed';src=OUT/'prospective_uncertainty_summary_full.csv'
d=pd.read_csv(src)
rows=[]
order=['exact_dose_reference','baseline','no_common_dose','no_relative_dose','no_oxygen','no_survival','improved_dose','improved_all','all_zero']
for sc in order:
    r=d[(d.scenario==sc)&(d.survival_rel_SE==.01)&(d.true_delta_percent==-2)].iloc[0]
    rows.append(dict(group='scenario',scenario=sc,survival_rel_SE=.01,true_delta_percent=-2,q2_5_percent=r.q2_5_percent,q97_5_percent=r.q97_5_percent,half_width_pp=r.half_width_pp))
for se in (.01,.03,.05):
    r=d[(d.scenario=='baseline')&(d.survival_rel_SE==se)&(d.true_delta_percent==-2)].iloc[0]
    rows.append(dict(group='survival_precision',scenario='baseline',survival_rel_SE=se,true_delta_percent=-2,q2_5_percent=r.q2_5_percent,q97_5_percent=r.q97_5_percent,half_width_pp=r.half_width_pp))
for td in (-2,-5,-10):
    r=d[(d.scenario=='baseline')&(d.survival_rel_SE==.01)&(d.true_delta_percent==td)].iloc[0]
    rows.append(dict(group='truth_sensitivity',scenario='baseline',survival_rel_SE=.01,true_delta_percent=td,q2_5_percent=r.q2_5_percent,q97_5_percent=r.q97_5_percent,half_width_pp=r.half_width_pp))
pd.DataFrame(rows).to_csv(OUT/'prospective_uncertainty_summary_manuscript.csv',index=False)
