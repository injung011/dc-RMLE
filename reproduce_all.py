"""Reproduce manuscript numerical outputs and figures.

Default mode is quick: deterministic analyses are recomputed and archived final
Monte Carlo summaries are used for the figures. Use --full to rerun the N=1000
retrospective and prospective Monte Carlo analyses, including the lipid-perturbation
recovery diagnostic, before regenerating all outputs.
"""
import argparse, os, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(rel, env=None, args=None):
    cmd=[sys.executable,str(ROOT/rel)] + (args or [])
    print('>', ' '.join(cmd)); subprocess.run(cmd,check=True,cwd=ROOT,env=env)

ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true',help='rerun N=1000 retrospective/prospective Monte Carlo analyses and lipid diagnostic');ap.add_argument('--workers',type=int,default=6);a=ap.parse_args()
if a.full:
    run('analyses/01_retrospective_fit_and_perturbation.py',args=['--n','1000','--beta-upper','10','--workers',str(a.workers)])
    run('analyses/01_retrospective_fit_and_perturbation.py',args=['--n','1000','--beta-upper','5','--workers',str(a.workers)])
    run('analyses/01b_prepare_beta_bound_summary.py')
    env=os.environ.copy();env['DC_RMLE_N_MC']='1000';env['DC_RMLE_WORKERS']=str(a.workers);run('analyses/04_prospective_uncertainty.py',env=env);run('analyses/04b_prepare_manuscript_summary.py');run('analyses/04c_prospective_lipid_recovery.py',env=env)
run('analyses/02_condition_deletion.py');run('analyses/03_oxygen_sensitivity.py');run('analyses/05_pratx_kapp_transfer.py')
run('analyses/06_export_manuscript_tables.py')
for f in ['make_figure1.py','make_figure2.py','make_figure3.py','make_figure4.py','make_figureS1.py']: run(Path('figures')/f)
print('Done. Figures: results/figures/ ; numerical outputs: data/processed/')
