import subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for f in ['make_figure1.py','make_figure2.py','make_figure3.py','make_figure4.py','make_figureS1.py']:
    subprocess.run([sys.executable,str(ROOT/'figures'/f)],cwd=ROOT,check=True)
print('Figures reproduced in results/figures/')
