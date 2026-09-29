"""Export manuscript-facing convenience tables from canonical processed outputs."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parents[1]; PROC=ROOT/'data'/'processed'; OUT=ROOT/'results'/'tables'; OUT.mkdir(parents=True,exist_ok=True)
mapping={
'retrospective_reported_only_summary.csv':'Table2_retrospective_summary.csv',
'Table3_ridge_points.csv':'Table3_hypCONV_deletion_ridge.csv',
'LQ_benchmark.csv':'TableS4_LQ_benchmark.csv',
'oxygen_trajectory_sensitivity.csv':'TableS5_oxygen_trajectory_sensitivity.csv',
'beta_bound_sensitivity.csv':'Supplementary_S3_4_beta_bound_sensitivity.csv',
'prospective_uncertainty_summary_manuscript.csv':'Prospective_uncertainty_summary.csv',
'TableS11_lipid_perturbation_recovery.csv':'TableS11_lipid_perturbation_recovery.csv',
'PratxKapp_transfer_summary.csv':'Table5_PratxKapp_transfer.csv',
'PratxKapp_alpha0_R_compensation.csv':'TableS9_PratxKapp_alpha0_compensation.csv',
'PratxKapp_consistent_Dmax_scan.csv':'TableS10_PratxKapp_Dmax_scan.csv',
'Figure2C_RMSE_minimum.csv':'Figure2C_RMSE_minimum.csv'}
missing=[]
for a,b in mapping.items():
    src=PROC/a
    if not src.exists(): missing.append(str(src.relative_to(ROOT)))
    else: shutil.copy2(src,OUT/b)
if missing: raise FileNotFoundError('Missing canonical outputs: '+', '.join(missing))
print(f'Exported {len(mapping)} manuscript-facing tables to results/tables/')
