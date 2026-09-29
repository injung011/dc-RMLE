# Supplementary insert: recovery of the imposed lipid perturbation

## Suggested subsection

### Recovery of the imposed lipid-directed perturbation

As an internal check that the synthetic lipid-directed perturbation generated a recoverable experimental contrast, the four independently fitted lipid-associated nuisance parameters were examined across the prospective Monte Carlo ensembles. The synthetic generator set each suppressed-state lipid contribution to 30% of its corresponding vehicle value, but this ratio was not supplied to the inverse fit. To avoid encoding the generating suppression ratio through the optimizer start point, all four lipid-associated terms were initialized to the same value (0.10) in this diagnostic.

Under the exact-dose reference, the median recovered suppressed-to-vehicle ratios were 0.258 for the CONV lipid term and 0.206 for the FLASH lipid term; the corresponding 2.5–97.5 percentile ranges were 0.087–0.423 and approximately 0–0.458. Under the baseline uncertainty scenario, the median ratios were 0.263 and 0.220, with ranges of 0.082–0.427 and approximately 0–0.460, respectively. In all 1000 realizations of both scenarios, each fitted suppressed-state lipid term was smaller than its corresponding vehicle term. In the zero-error implementation check, the generating ratio of 0.300 was recovered exactly. Thus, the design robustly recovered the direction of the imposed lipid perturbation, whereas the exact 70% suppression magnitude remained only weakly determined.

## Suggested Table S11 caption

**Table S11. Recovery of the independently fitted lipid-associated nuisance terms in the prospective synthetic experiment.** The generator set the suppressed-state lipid contribution to 30% of the corresponding vehicle value for both CONV and FLASH, but this ratio was not supplied to the inverse fit. Values are medians with 2.5–97.5 percentile ranges over 1000 Monte Carlo realizations. The final column gives the fraction of realizations in which the fitted suppressed-state term was smaller than the corresponding vehicle term.

| Scenario | Mode | Vehicle Γ, median (2.5–97.5%) | Suppressed Γ, median (2.5–97.5%) | Suppressed/vehicle ratio, median (2.5–97.5%) | Fraction with Γsupp < Γvehicle |
|---|---|---|---|---|---:|
| Exact-dose reference | CONV | 0.124 (0.097–0.163) | 0.0317 (0.0086–0.0687) | 0.258 (0.087–0.423) | 1.000 |
| Exact-dose reference | FLASH | 0.0905 (0.0683–0.1349) | 0.0188 (~0–0.0618) | 0.206 (~0–0.458) | 1.000 |
| Baseline uncertainty | CONV | 0.124 (0.096–0.164) | 0.0328 (0.0081–0.0694) | 0.263 (0.082–0.427) | 1.000 |
| Baseline uncertainty | FLASH | 0.0919 (0.0675–0.1349) | 0.0201 (~0–0.0622) | 0.220 (~0–0.460) | 1.000 |

Generator values: ΓC,vehicle = 0.130429, ΓF,vehicle = 0.102507, ΓC,suppressed = 0.039129, and ΓF,suppressed = 0.030752; the true suppressed/vehicle ratio is 0.300 for both modes.
