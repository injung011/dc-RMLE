# dc-RMLE identifiability and parameter-recovery analyses

Reproducibility package for the manuscript **“What can be inferred from FLASH clonogenic survival? Parameter recovery and the role of irradiation conditions in dc-RMLE.”**

## Scope

This repository reproduces the numerical analyses in the manuscript:

- joint reference fitting of four reconstructed DU145 survival curves and the reported-only reconstruction-perturbation ensemble;
- independent non-negative LQ descriptive fits used in Figure 1;
- one-condition deletion diagnostics, including the exact hypoxia-CONV reduced-design rescaling ridge;
- sensitivity to the prescribed FLASH oxygen-depletion trajectory;
- prospective direct-oxygen + lipid-perturbation Monte Carlo recovery analysis with explicit measurement uncertainty;
- diagnostic recovery of the independently fitted vehicle and lipid-suppressed nuisance terms, without supplying the generating 30% suppressed/vehicle ratio to the inverse fit;
- Pratx–Kapp inverse-problem transfer analysis, including the depletion-breakpoint dose scan.

## Interpretation

The retrospective percentile ranges are **reconstruction-perturbation sensitivity ranges**, not formal confidence intervals. The hypoxia-CONV deletion ridge is an exact algebraic rescaling of the reduced design, not an optimizer artifact. The prospective study is a synthetic design study under explicitly stated measurement assumptions. The Pratx–Kapp calculation is an inverse-problem transfer check and does not reassess the original forward-model study.

## Data sources

`data/raw/DU145_reconstructed.csv` contains 26 non-zero-dose clonogenic-survival points reconstructed from the published DU145 curves of Adrian et al. (2020), together with graphically reconstructed error-bar values where they could be resolved. These are **digitized/reconstructed values, not the original replicate-level data**. Missing graphical errors are retained as missing and are not perturbed in the reported-only reconstruction ensemble. Users of these reconstructed data should cite the source publication:

> Adrian G, Konradsson E, Lempart M, Bäck S, Ceberg C and Petersson K 2020. *The FLASH effect depends on oxygen concentration.* Br. J. Radiol. 93 20190702.

The Pratx–Kapp transfer analysis implements the oxygen-depletion/OER functional form of Pratx and Kapp (2019). The values `R = 1.63`, `phi = 0.26 mmHg^-1`, and `L_ROD = 0.42 mmHg Gy^-1` follow that formulation; `alpha0 = 0.15 Gy^-1` is a numerical reference value used for the inverse-problem check in this repository rather than a fitted parameter reported by the original study. Please cite:

> Pratx G and Kapp D S 2019. *A computational model of radiolytic oxygen depletion during FLASH irradiation and its effect on the oxygen enhancement ratio.* Phys. Med. Biol. 64 185005.

## Environment

The final package was checked under Python 3.13 with NumPy 2.3, pandas 2.2, SciPy 1.17, Matplotlib 3.10 and Pillow 12. Install compatible versions with:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Reproduce outputs

For a fast reconstruction of all deterministic analyses and all manuscript figures using the archived final Monte Carlo summaries:

```bash
python reproduce_all.py
```

To rerun the `N = 1000` retrospective and prospective Monte Carlo analyses, including the lipid-perturbation recovery diagnostic, and then regenerate all figures:

```bash
python reproduce_all.py --full --workers 6
```

The full run is substantially slower because it performs repeated bounded nonlinear least-squares fits.

To regenerate figures only:

```bash
python reproduce_figures.py
```

## Random seeds

- retrospective reported-only reconstruction perturbation: `20260921`
- prospective uncertainty propagation and lipid-perturbation recovery diagnostic: `20260920`
- prospective scenario comparisons use common random numbers across matched uncertainty scenarios;
- the deletion contour and exact rescaling ridge are deterministic and do not use a random seed.

## Data and outputs

`data/processed/` contains the canonical numerical outputs used by the analyses and figures. The distributed package also contains the final Monte Carlo sample files used to compute the archived summaries; they can be regenerated with `--full`.

`results/tables/` contains manuscript-facing convenience exports synchronized from the canonical processed outputs by `analyses/06_export_manuscript_tables.py`. `results/figures/` contains PNG and vector PDF outputs.

See `docs/MANUSCRIPT_ITEM_MAP.md` for the direct mapping between manuscript figures/tables and their generating scripts.

## Key manuscript checks

A successful reproduction should recover, to rounding, the following checks.

### Retrospective reference fit and reconstruction perturbation

- reference dc-RMLE SSE ≈ `0.469806` and residual factor `Delta Ps/Ps ≈ -2.18%`;
- reported-only retrospective `Delta Ps/Ps` 2.5–97.5 percentile range ≈ `-17.96% to +5.66%`.

### Condition deletion and oxygen-trajectory sensitivity

- hypoxia-CONV deletion training SSE ≈ `0.294177` along the exact ridge;
- representative ridge scan: `R_Lambda ≈ 1.145–2.279` and `Delta Ps/Ps ≈ -12.9% to +73.4%`;
- invariant combinations along the exact reduced-design ridge are approximately
  - `qFLASH/(1 + GammaC) = 0.883312`,
  - `alpha2*(1 + GammaC) = 0.185616`,
  - `beta*(1 + GammaC)^2 = 1.100863`;
- oxygen-trajectory scan at `g0 = 0.10, 0.20, 0.40 mmHg Gy^-1`: `Delta Ps/Ps ≈ -8.50%, -6.19%, -2.18%`, respectively.

### Prospective uncertainty analysis

- exact-dose reference half-width ≈ `1.69 pp`;
- baseline dose-uncertainty half-width ≈ `4.97 pp`;
- removing the CONV–FLASH relative dose error reduces the half-width to ≈ `1.69 pp`;
- improved dose delivery gives ≈ `2.51 pp`;
- improved dose delivery plus improved oxygen precision gives ≈ `2.31 pp`;
- with all measurement errors set to zero, the generating value is recovered exactly: `Delta Ps/Ps = -2.00%` with zero width.

### Prospective lipid-perturbation recovery diagnostic

The synthetic generator sets each suppressed-state lipid contribution to `0.30` times its corresponding vehicle value, but this ratio is not constrained in the inverse fit. The four lipid nuisance parameters are fitted independently. For this diagnostic, all four lipid terms are initialized to the same neutral value (`0.10`) so the generating suppression ratio is not encoded in the optimizer start point.

- exact-dose reference median suppressed/vehicle ratios: `0.258` (CONV) and `0.206` (FLASH);
- baseline-uncertainty median suppressed/vehicle ratios: `0.263` (CONV) and `0.220` (FLASH);
- in all `1000/1000` realizations of both scenarios, the fitted suppressed-state lipid term was smaller than its corresponding vehicle term for both CONV and FLASH;
- the zero-error diagnostic recovers the generating ratio exactly (`0.300`).

These checks support recovery of the **direction** of the imposed lipid perturbation; the broad ratio intervals show that the exact 70% suppression magnitude is not tightly recovered.

### Pratx–Kapp transfer analysis

- fixed-`alpha0` compensation: `alpha0(1 + R) = 0.394500` is conserved for `alpha0 = 0.10–0.39 Gy^-1`;
- **main-text Table 5 designs**: scaled Jacobian condition number ≈ `571.5` (CONV only), `584.2` (FLASH to 8 Gy using four dose points: 2, 4, 6, 8 Gy), `37.9` (FLASH to 20 Gy), and `9.67` (with a separate anoxic reference);
- **Supplementary Table S10 scan**, using eight equally spaced FLASH doses from 2 Gy to `D_max`: ≈ `595.8` at 8 Gy, `600.1` at the 9.52 Gy breakpoint, `491.2` at 10 Gy, and `36.7` at 20 Gy.

The Table 5 and Table S10 values at 8 Gy differ because the two calculations use different FLASH dose grids.

## Release-integrity note

`SHA256SUMS.txt` verifies the distributed release files. Regenerating figures can change byte-level PDF metadata even when the plotted numerical content is unchanged; numerical checks should therefore use the CSV outputs and the rounded values listed above.

## How to cite

If you use this code or the reconstructed data, please cite the manuscript and this repository:

- Manuscript: `[citation to be added on publication]`
- Repository archive: `[DOI to be added]`

## License

The software and analysis code in this repository are released under the MIT License. See the `LICENSE` file for details.

Unless otherwise stated, this license applies to the software and analysis code provided in this repository and does not alter the copyright or licensing terms of third-party data, publications, or other externally sourced materials.
