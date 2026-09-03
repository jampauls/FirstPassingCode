# Test Scripts

This directory contains MATLAB R2020b regression tests, verification
scripts, and OWNS analysis studies for the first-transition codebase.
Most scripts add the project `src` and `tests` directories automatically.
Run them from MATLAB with the project root on the path, for example:

```matlab
cd('/home/jampauls/Documents/trimmed_exactFirstPassingCode')
run('tests/runSyntheticRegressionTests.m')
```

The OWNS scripts require the OWNS solution MAT-file configured in their
`cfg.owns.dataFile` fields. The checked-in default is an absolute path under
`/data2/jampauls/...`; update that setting when the data live elsewhere.
Several OWNS studies save MAT reports to the project `results` directory and
open diagnostic figures.

## Regression and Verification

| Script | Purpose | Inputs and checks |
| --- | --- | --- |
| `runSyntheticRegressionTests.m` | End-to-end synthetic first-transition regression suite. | Builds the rotating 2D problem; verifies analytical matrices, circular quadrature convergence, trajectory Monte Carlo agreement, local-event consistency, probability conservation, and a terminal CDF regression target. |
| `runHermiteMaximumTests.m` | Verifies continuous cubic-Hermite maximum detection. | Checks exact cubic maximization, exact and finite-difference derivatives, agreement with grid and dense-grid CDFs, then plots diagnostics. |
| `runRQMCTests.m` | Verifies scrambled-Sobol angular integration. | Compares RQMC resolutions against a high-resolution circle rule; checks statistical agreement, local-event consistency, repeatability, and independent replicates. |
| `runBivariateQuadraticProbabilityTests.m` | Tests 2D quadratic-form inversion and facet exclusion. | Covers independent forms, identical facets, proportional ordered forms, and impossible ordered exclusion. The independent-form assertion is intentionally commented out; the remaining checks are active. |

## OWNS End-to-End Studies

| Script | Purpose | Main output |
| --- | --- | --- |
| `runOWNSSmokeTest.m` | Small end-to-end OWNS calculation. | Loads/compresses OWNS data, builds derivatives, computes an RQMC CDF, validates it, and plots energy, CDF, survival, and memory correction. |
| `runOWNSReductionStudy.m` | Uncoupled transition-oriented rank convergence. | Compares integrated-basis ranks using RQMC and plots terminal CDFs, rank changes, and matrix residuals. |
| `runOWNSCoupledReductionStudy.m` | Coupled rank convergence with common RQMC directions. | Computes paired rank differences, matrix residuals, and an automatic rank recommendation. |
| `runOWNSThresholdStudy.m` | Threshold and inlet-amplitude sensitivity. | Reuses an RQMC energy envelope to evaluate transition CDFs over threshold and amplitude grids. |

## OWNS Deterministic Envelope Studies

| Script | Purpose | Main output |
| --- | --- | --- |
| `runOWNSDeterministicLoewnerBoundStudy.m` | Computes deterministic first-transition brackets using a positive-variation Loewner envelope. | Saves `results/owns_deterministic_loewner_bound_study.mat` and reports scalar inversion and envelope diagnostics. |
| `runOWNSDeficitLoewnerEnvelopeStudy.m` | Compares positive-variation and sequential deficit Loewner envelopes. | Saves `results/owns_deficit_loewner_envelope_study.mat`; validates dominance and monotonicity before evaluating bounds. |
| `runOWNSPositiveIncrementStudy.m` | Diagnoses whether growth lies in a low-dimensional stochastic subspace. | Analyzes positive increments, complement growth, and optional directional record-growth sampling. |

## OWNS Matrix-Feature Studies

| Script | Purpose | Main output |
| --- | --- | --- |
| `runOWNSMatrixFeatureStudy.m` | Measures low-rank structure in matrix trajectories. | Analyzes linear, affine, and positive-increment families; saves `results/owns_matrix_feature_study.mat`. |
| `runOWNSMatrixFeatureProbabilityStudy.m` | Validates feature truncations against the full first-transition CDF. | Uses coupled RQMC validation and saves `results/owns_matrix_feature_probability_study.mat`. It reuses the matrix-feature report when available. |

## Notes

- Synthetic tests are the quickest checks and do not depend on the OWNS data file.
- The RQMC and dense-grid checks can be computationally expensive because they use large direction counts, multiple replicates, or dense station grids.
- OWNS scripts assume a constant threshold where required and stop with an error when that assumption is violated.