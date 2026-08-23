# OWNS First-Transition Probability

MATLAB implementation of a memory-aware first-transition method for Gaussian disturbances propagated by the One-Way Navier–Stokes (OWNS) equations.

The code computes the probability that a disturbance-energy threshold has been crossed at or before each streamwise station. Unlike a local energy-exceedance probability, the first-transition probability retains the upstream history of each stochastic realization.

## Mathematical overview

The propagated disturbance is represented in a common set of inlet stochastic coordinates:

$$
q(x)=B(x)w, \quad w \sim N(0, I_r).
$$

The default implementation uses **real Gaussian coefficients**, even though the Fourier-transformed OWNS states and propagated basis vectors can be complex.

For an energy weight matrix $H(x)$,

$$
e(x)=q(x)^*H(x)q(x)=w^T A(x)w,
$$

where

$$
A(x)=\mathrm{Re}[B(x)^*H(x)B(x)].
$$

Given a positive transition threshold $e_{\mathrm{thres}}(x)$, define

$$
G(x)=\frac{A(x)}{e_{\mathrm{thres}}(x)}.
$$

The first-transition location is

$$
X_{\mathrm{tr}}=\inf\{x: w^T G(x) w \ge 1\}.
$$

Using the decomposition

$$
w=Ru, \quad R^2 \sim \chi_r^2, \quad u \sim \mathrm{Unif}(S^{r-1}),
$$

define the directional running maximum

$$
m_u(x)=\max_{\xi \le x} u^T G(\xi)u.
$$

The first-transition cumulative distribution function is then

$$
F_{X_{\mathrm{tr}}}(x)=E_u\left[\bar F_{\chi_r^2}\left(\frac{1}{m_u(x)}\right)\right].
$$

The running maximum provides streamwise memory: once a realization crosses the threshold, it remains classified as transitioned downstream.

## Current capabilities

The code currently supports:

- synthetic rotating quadratic-form test problems;
- precomputed OWNS solution files;
- real Gaussian stochastic coefficients;
- proper-complex Gaussian coefficients as a sensitivity model;
- complex Fourier-state factors;
- pressure-excluded OWNS energy norms;
- grid-based and cubic-Hermite streamwise maximum detection;
- deterministic circular quadrature for two stochastic dimensions;
- random angular integration;
- scrambled Sobol’ randomized quasi-Monte Carlo;
- replicate-based RQMC standard-error estimates;
- transition-oriented stochastic dimension reduction;
- integrated, eigenspace-union, hybrid, and greedy bases;
- coupled RQMC rank-convergence studies;
- local energy-exceedance probabilities;
- independent characteristic-function inversion for local quadratic forms;
- survival probabilities;
- first-transition CDFs;
- interval transition probabilities;
- right-censoring probabilities;
- transition quantiles;
- local-versus-history-aware memory corrections;
- threshold and inlet-amplitude sensitivity studies;
- interactive OWNS file selection with `uigetfile`;
- real-versus-proper-complex Gaussian convention studies.

## MATLAB requirements

The target version is:

- MATLAB R2020b

The following toolbox is required for scrambled Sobol’ RQMC:

- Statistics and Machine Learning Toolbox

In particular, the RQMC implementation uses:

```matlab
sobolset
scramble
```

Most chi-square and gamma probabilities are evaluated through MATLAB’s built-in `gammainc`, avoiding a dependency on `chi2cdf`.

A graphical MATLAB session is required only when using:

```matlab
cfg.dataSource = 'ownsUserSelect';
```

For batch or headless execution, use:

```matlab
cfg.dataSource = 'ownsMatFile';
```

## Repository structure

```text
project_root/
├── README.md
├── run_first_transition.m
│
├── src/
│   ├── angular/
│   │   ├── computeProperComplexRQMCCDF.m
│   │   ├── computeRQMCCDF.m
│   │   ├── computeRQMCCDFConverged.m
│   │   ├── evaluateAngularReplicate.m
│   │   ├── generateAngularRule.m
│   │   ├── generateProperComplexAngularRule.m
│   │   ├── generateScrambledSobolNormals.m
│   │   ├── nestedNormalsToDirections.m
│   │   └── validateAngularRule.m
│   │
│   ├── data/
│   │   ├── convertFactorOrientation.m
│   │   ├── extractOWNSCoordinate.m
│   │   ├── inspectOWNSData.m
│   │   ├── loadOWNSProblem.m
│   │   ├── loadProblemData.m
│   │   ├── makeSyntheticProblem.m
│   │   ├── runOWNSWorkflow.m
│   │   ├── selectOWNSDataFile.m
│   │   └── validateProblemData.m
│   │
│   ├── diagnostics/
│   │   └── computeDirectionalDiagnostics.m
│   │
│   ├── matrices/
│   │   ├── buildAFromOWNS.m
│   │   ├── buildGFamily.m
│   │   ├── buildGFromAFamily.m
│   │   ├── buildGprimeFamily.m
│   │   ├── validateGFamily.m
│   │   └── validateGprimeFamily.m
│   │
│   ├── maxima/
│   │   ├── computeRunningMaxAdaptive.m
│   │   ├── computeRunningMaxGrid.m
│   │   ├── computeRunningMaxHermite.m
│   │   ├── cubicHermiteIntervalMaximum.m
│   │   └── validateRunningMax.m
│   │
│   ├── plotting/
│   │   └── plotTransitionResults.m
│   │
│   ├── probability/
│   │   ├── buildRQMCEnergyEnvelope.m
│   │   ├── buildThresholdSweepValues.m
│   │   ├── buildTransitionThreshold.m
│   │   ├── computeLocalExceedance.m
│   │   ├── computeLocalExceedanceCF.m
│   │   ├── computeLocalExceedanceFromGains.m
│   │   ├── computeProperComplexLocalExceedanceFromGains.m
│   │   ├── computeProperComplexTrajectoryMonteCarlo.m
│   │   ├── computeProperComplexTransitionCDF.m
│   │   ├── computeTrajectoryMonteCarlo.m
│   │   ├── computeTransitionCDF.m
│   │   ├── computeTransitionQuantiles.m
│   │   ├── evaluateEnergyEnvelope.m
│   │   ├── evaluateThresholdAmplitudeFamily.m
│   │   ├── generalizedQuadraticFormTail.m
│   │   ├── solveQuadraticFormThreshold.m
│   │   ├── solveThresholdForTargetProbability.m
│   │   └── validateProbabilityCurves.m
│   │
│   └── reduction/
│       ├── applyComplexReductionToMatrixFamily.m
│       ├── applyReductionToFactors.m
│       ├── applyReductionToMatrixFamily.m
│       ├── buildComplexIntegratedBasis.m
│       ├── buildEigenspaceSnapshotBasis.m
│       ├── buildIntegratedTransitionMatrix.m
│       ├── buildTransitionBasis.m
│       ├── computeCoupledComplexRankStudy.m
│       ├── computeCoupledRankStudy.m
│       ├── computeReductionDiagnostics.m
│       ├── computeTrapezoidalWeights.m
│       ├── enrichTransitionBasisGreedy.m
│       ├── selectRankFromCoupledStudy.m
│       ├── selectTransitionStations.m
│       ├── validateComplexReductionBasis.m
│       └── validateReductionBasis.m
│
├── examples/
│   └── plotSyntheticProblem.m
│
├── tests/
│   ├── runHermiteMaximumTests.m
│   ├── runOWNSCoupledReductionStudy.m
│   ├── runOWNSGaussianConventionStudy.m
│   ├── runOWNSLocalProbabilityCheck.m
│   ├── runOWNSProperComplexRankStudy.m
│   ├── runOWNSReductionStudy.m
│   ├── runOWNSSmokeTest.m
│   ├── runOWNSThresholdStudy.m
│   ├── runProperComplexRegressionTests.m
│   ├── runQuadraticFormTailTests.m
│   ├── runRQMCTests.m
│   └── runSyntheticRegressionTests.m
│
└── data/
    └── .gitkeep
```

Some files shown above may remain stubs until their corresponding capability is completed. In particular, the adaptive maximum finder and direct OWNS marching interface are not part of the current production path.

## Getting started

From the project root, launch MATLAB and run:

```matlab
run_first_transition
```

The script automatically adds the source, example, and data directories to the MATLAB path.

For tests, add the test directory if it is not already included:

```matlab
addpath(genpath(fullfile(pwd, 'tests')));
```

## Synthetic example

The synthetic example is the recommended first check after cloning the repository.

Use:

```matlab
cfg.dataSource = 'synthetic';
cfg.gaussianConvention = 'real';

cfg.angular.method = 'circle';
cfg.angular.numDirections = 4096;

cfg.maxDetection.method = 'hermite';
cfg.reduction.method = 'none';
```

Then run:

```matlab
run_first_transition
```

For the current two-dimensional synthetic problem, the expected Hermite-corrected terminal probability is approximately:

```text
F(xmax) = 3.461875e-01
```

The exact reported digits depend on angular resolution.

## OWNS input format

The current OWNS adapter expects a MAT-file containing a structure named:

```matlab
solution
```

with fields compatible with:

```matlab
solution.q_in
solution.q
solution.gram_W
solution.x
solution.y
solution.z
solution.eta
solution.xi
solution.beta
solution.w
solution.fp
```

For the current data:

```matlab
Ny = 201;
Nx = 268;
r  = 141;
```

The propagated factor has size:

```matlab
size(solution.q)
% 6*Ny by Nx by r
```

The state uses variable-major ordering:

```text
[rho(1:Ny),
 u(1:Ny),
 v(1:Ny),
 w(1:Ny),
 T(1:Ny),
 p(1:Ny)]
```

Pressure is excluded from the disturbance-energy norm. The energy-bearing state is therefore:

```matlab
energyIndices = 1:(5*Ny);
```

At station `n`, the synthesis factor is:

```matlab
Bn = squeeze(solution.q(energyIndices, n, :));
```

The diagonal energy weights are:

```matlab
weights = full(solution.gram_W(:, n));
```

The stochastic-space energy matrix is:

```matlab
weightedB = bsxfun(@times, sqrt(weights), Bn);
Acomplex = weightedB' * weightedB;
Areal = real(Acomplex);
Areal = 0.5 * (Areal + Areal.');
```

The OWNS columns are assumed to remain in one-to-one correspondence with the inlet columns after the initial downstream-mode projection.

The temporal frequency is read only from:

```matlab
solution.w
```

The field:

```matlab
solution.qCart.w
```

denotes spanwise velocity and must not be interpreted as temporal frequency.

## Selecting an OWNS data file

### Fixed path

```matlab
cfg.dataSource = 'ownsMatFile';

cfg.owns.dataFile = ...
    '/path/to/owns_solution.mat';
```

### Interactive selection

```matlab
cfg.dataSource = 'ownsUserSelect';

cfg.owns.userSelectDirectory = ...
    '/directory/containing/owns/files';

cfg.owns.userSelectTitle = ...
    'Select an OWNS solution MAT-file';

cfg.owns.userSelectDefaultFile = '';
```

The selected path is stored in:

```matlab
cfg.owns.dataFile
meta.selectedDataFile
```

## Streamwise coordinate

The recommended first-passage coordinate is wall arc length:

```matlab
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;
```

Available coordinate options include:

```text
wallArcLength
referenceArcLength
wallX
referenceX
xi
```

For curved geometries, wall arc length is computed from the physical `x`, `y`, and optionally `z` coordinates.

## Transition threshold

A validated physical threshold is not yet available. The code supports provisional data-derived thresholds and user-specified values.

### Automatic development threshold

```matlab
cfg.threshold.method = 'auto';
cfg.threshold.autoMaxMeanFactor = 1;
```

This uses:

\[
e_{\mathrm{thres}}
=
\max_x \mathbb E[e(x)].
\]

This is intended only for numerical development.

### Specified constant threshold

```matlab
cfg.threshold.method = 'specifiedScalar';
cfg.threshold.value = 2.5e-4;
```

### Specified streamwise threshold

```matlab
cfg.threshold.method = 'specifiedVector';
cfg.threshold.vector = thresholdVector;
```

### Inlet exceedance quantile

```matlab
cfg.threshold.method = 'initialQuantileCF';
cfg.threshold.exceedanceProbability = 1e-3;
```

This uses deterministic characteristic-function inversion of the inlet quadratic-form distribution.

Any automatically generated threshold should be treated as provisional until calibrated or justified using DNS or another physical transition criterion.

## Stochastic coefficient convention

### Default: real Gaussian coefficients

The primary model is:

```matlab
w = randn(r, 1);
qRealization = B * w;
physicalField = real(qRealization);
```

Use:

```matlab
cfg.gaussianConvention = 'real';
```

The radial law is:

\[
R^2\sim\chi_r^2.
\]

### Sensitivity model: proper-complex coefficients

The optional sensitivity model is:

```matlab
z = (randn(r,1) + 1i*randn(r,1)) / sqrt(2);
qRealization = B * z;
```

Use the proper-complex probability functions and complex-specific reduction basis.

The radial law is:

\[
R^2\sim\operatorname{Gamma}(r,1).
\]

The proper-complex model represents additional independent random phase variation in each latent coefficient. It is not currently the primary physical model.

## Streamwise maximum detection

### Stored-grid maximum

```matlab
cfg.maxDetection.method = 'grid';
```

This checks only stored streamwise stations and generally provides a lower approximation to continuous first-transition probability.

### Cubic-Hermite maximum

```matlab
cfg.maxDetection.method = 'hermite';
```

This uses endpoint values and derivatives of each directional gain to locate interior maxima between stations.

For precomputed OWNS data, derivatives are currently estimated with a nonuniform three-point polynomial stencil.

### Adaptive method

```matlab
cfg.maxDetection.method = 'adaptive';
```

This branch is reserved for future work and is not currently part of the verified production workflow.

## Angular integration

### Two-dimensional deterministic rule

For \(r=2\):

```matlab
cfg.angular.method = 'circle';
cfg.angular.numDirections = 4096;
```

### Scrambled Sobol’ RQMC

For moderate or high stochastic dimension:

```matlab
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 4096;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 1;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;
```

RQMC replicates are evaluated using common local and first-transition directional samples, preserving:

\[
p_{\mathrm{local}}(x)
\leq
F_{X_{\mathrm{tr}}}(x)
\]

replicate by replicate.

## Transition-oriented dimension reduction

The full OWNS stochastic rank is currently \(r=141\). A fixed reduced basis is used over the complete streamwise domain.

### Integrated basis

```matlab
cfg.reduction.method = 'integrated';
cfg.reduction.targetRank = 30;
cfg.reduction.integrationWeight = 'uniform';
```

The integrated matrix is:

\[
K
=
\int G(x)\,dx.
\]

Its leading eigenvectors define the reduced stochastic basis.

### Other available basis methods

```text
eigUnion
hybrid
greedy
givenBasis
none
```

For the current real-coefficient OWNS dataset, a coupled RQMC rank study indicated that rank 30 satisfies an approximate absolute complete-CDF tolerance of \(5\times10^{-3}\) under the provisional threshold.

Rank adequacy should be checked again if:

- the dataset changes materially;
- the threshold becomes streamwise varying;
- the energy norm changes;
- or the Gaussian coefficient convention changes.

## Local energy probabilities

Two local probability methods are available.

### Same-angular-rule estimate

This uses the same RQMC directions as the first-transition calculation:

```matlab
cfg.localProbability.method = 'sameAngular';
```

### Independent characteristic-function inversion

```matlab
cfg.localProbability.method = 'characteristicFunction';
```

or:

```matlab
cfg.localProbability.method = 'both';
```

The deterministic method evaluates the local generalized quadratic-form tail using characteristic-function inversion and is independent of angular RQMC.

## Threshold and amplitude studies

For a constant threshold and inlet amplitude multiplier \(\epsilon\),

\[
F_{X_{\mathrm{tr}}}(x)
=
\mathbb E_u
\left[
    \overline F_{\chi_r^2}
    \left(
        \frac{e_{\mathrm{thres}}}
        {\epsilon^2 M_u(x)}
    \right)
\right].
\]

A reusable unnormalized energy envelope permits many thresholds and amplitudes to be evaluated without repeating the directional maximum calculation.

Enable the main-script sweep with:

```matlab
cfg.thresholdSweep.enable = true;
```

Representative options are:

```matlab
cfg.thresholdSweep.factors = logspace(-0.5, 0.5, 9);
cfg.thresholdSweep.amplitudes = [0.5, 0.75, 1.0, 1.25, 1.5];
cfg.thresholdSweep.referenceMethod = 'maximumMeanEnergy';
```

A target terminal probability can also be used to solve for a constant threshold.

## Tests

Run the following tests from the MATLAB command window.

### Synthetic regression

```matlab
runSyntheticRegressionTests
```

Verifies:

- construction of the synthetic \(G(x)\);
- angular convergence;
- trajectory Monte Carlo agreement;
- local-versus-first-transition consistency;
- probability conservation.

### Hermite maximum detection

```matlab
runHermiteMaximumTests
```

Verifies:

- exact reproduction of cubic maxima;
- exact synthetic derivatives;
- Hermite versus stored-grid maxima;
- agreement with a dense streamwise reference;
- finite-difference derivative fallback.

### RQMC integration

```matlab
runRQMCTests
```

Verifies:

- scrambled Sobol’ convergence;
- deterministic reproducibility under fixed seeds;
- agreement with a deterministic circular reference;
- replicate uncertainty.

### OWNS smoke test

```matlab
runOWNSSmokeTest
```

Verifies the complete OWNS loading, energy compression, threshold, Hermite, and RQMC path.

### Real-model rank convergence

```matlab
runOWNSCoupledReductionStudy
```

Uses common RQMC coordinates across nested ranks to estimate paired reduction error.

### Proper-complex verification

```matlab
runProperComplexRegressionTests
```

Verifies:

- analytical scalar proper-complex probabilities;
- proper-complex energy moments;
- direct trajectory Monte Carlo agreement.

### Proper-complex rank convergence

```matlab
runOWNSProperComplexRankStudy
```

Checks complex-specific transition-oriented bases and full-rank convergence.

### Generalized quadratic-form tails

```matlab
runQuadraticFormTailTests
```

Verifies deterministic local tail probabilities and threshold inversion.

### OWNS local probability comparison

```matlab
runOWNSLocalProbabilityCheck
```

Compares angular RQMC local probabilities against independent characteristic-function inversion.

### Threshold sensitivity

```matlab
runOWNSThresholdStudy
```

Evaluates constant-threshold and inlet-amplitude sensitivity using a reusable energy envelope.

## Representative verified results

### Synthetic two-dimensional test

Endpoint-only circular result:

```text
F(xmax) = 0.3461435512
```

Hermite-corrected result:

```text
F(xmax) = 0.3461875
```

### Current OWNS dataset with provisional threshold

Using the real-coefficient model and the provisional automatic threshold, the full-rank transition probability was approximately:

```text
F(xmax) ≈ 0.447
```

Using proper-complex coefficients produced a larger probability of approximately:

```text
F(xmax) ≈ 0.483
```

These values are for numerical development only because the physical transition threshold has not yet been calibrated.

## Current limitations

The following capabilities are not currently part of the verified production path:

- direct execution of the OWNS spatial march;
- nonzero distributed forcing;
- adaptive insertion of new physical OWNS stations;
- nonlinear disturbance propagation;
- DNS threshold calibration;
- broadband joint reconstruction across multiple \((\omega,\beta)\) pairs;
- directory-wide processing of all spectral pairs;
- physical-field plotting and reconstruction inside this repository;
- automatic parallel RQMC execution.

OWNS runs and physical disturbance visualization are expected to remain external to this codebase.

## Planned future work

Current priorities are:

1. finish verification of deterministic generalized quadratic-form tails;
2. add configurable higher-order derivative stencils for \(G'(x)\);
3. automate RQMC convergence for production calculations;
4. add optional parallel execution over RQMC replicates;
5. implement survivor-conditioned stochastic covariance;
6. implement transition-conditioned directional diagnostics;
7. design directory-wide processing of all \((\omega,\beta)\) pairs;
8. design joint reconstruction over the full spectral space.

The broadband reconstruction requires a deliberate modeling decision about spectral correlations, phase conventions, and whether transition is defined from summed energy or another joint spectral statistic.

## Reproducibility

For reproducible calculations, preserve:

- the complete configuration structure;
- selected OWNS file path;
- temporal frequency and spanwise wavenumber;
- stochastic rank and basis;
- threshold definition;
- maximum-detection method;
- RQMC point count;
- RQMC replicate count;
- all random seeds;
- MATLAB version;
- numerical tolerances.

The main result structure should include, at minimum:

```matlab
results.cfg
results.meta
results.prob
results.post
results.reductionInfo
```

## Data and repository policy

Large OWNS MAT-files should not be committed to Git. Store them externally and configure their paths at runtime.

Likewise, generated result files, figures, and temporary MATLAB files should generally be excluded unless intentionally added as small regression references.

## License

No license has yet been selected. Before publishing the repository, add a license file appropriate for the intended use and distribution policy.

## Citation

A formal citation has not yet been assigned. If this repository supports a publication or dissertation, add the corresponding citation and DOI here.

## Suggested .gitignore

```gitignore
# MATLAB temporary and autosave files
*.asv
*.m~
*.autosave

# MATLAB generated files
*.mex*
*.mlappinstall
*.mltbx

# Large data and result files
*.mat
*.h5
*.hdf5

# Generated output directories
results/
figures/
output/
tmp/
temp/

# Keep empty project data directory
!data/.gitkeep

# Profiling and coverage output
profile_results/
coverage/

# Operating-system files
.DS_Store
Thumbs.db

# Editor and IDE settings
.vscode/
.idea/

# Slurm or cluster output
slurm-*.out
*.o
*.e
```

If you intend to commit small MAT-files for regression tests, replace the broad `*.mat` rule with:

```gitignore
data/*.mat
results/*.mat
output/*.mat

!tests/reference_data/*.mat
```
referenceArcLength
wallX
referenceX
xi
```

For curved geometries, wall arc length is computed from the physical `x`, `y`, and optionally `z` coordinates.

## Transition threshold

A validated physical threshold is not yet available. The code supports provisional data-derived thresholds and user-specified values.

### Automatic development threshold

```matlab
cfg.threshold.method = 'auto';
cfg.threshold.autoMaxMeanFactor = 1;
```

This uses:

\[
e_{\mathrm{thres}}
=
\max_x \mathbb E[e(x)].
\]

This is intended only for numerical development.

### Specified constant threshold

```matlab
cfg.threshold.method = 'specifiedScalar';
cfg.threshold.value = 2.5e-4;
```

### Specified streamwise threshold

```matlab
cfg.threshold.method = 'specifiedVector';
cfg.threshold.vector = thresholdVector;
```

### Inlet exceedance quantile

```matlab
cfg.threshold.method = 'initialQuantileCF';
cfg.threshold.exceedanceProbability = 1e-3;
```

This uses deterministic characteristic-function inversion of the inlet quadratic-form distribution.

Any automatically generated threshold should be treated as provisional until calibrated or justified using DNS or another physical transition criterion.

## Stochastic coefficient convention

### Default: real Gaussian coefficients

The primary model is:

```matlab
w = randn(r, 1);
qRealization = B * w;
physicalField = real(qRealization);
```

Use:

```matlab
cfg.gaussianConvention = 'real';
```

The radial law is:

\[
R^2\sim\chi_r^2.
\]

### Sensitivity model: proper-complex coefficients

The optional sensitivity model is:

```matlab
z = (randn(r,1) + 1i*randn(r,1)) / sqrt(2);
qRealization = B * z;
```

Use the proper-complex probability functions and complex-specific reduction basis.

The radial law is:

\[
R^2\sim\operatorname{Gamma}(r,1).
\]

The proper-complex model represents additional independent random phase variation in each latent coefficient. It is not currently the primary physical model.

## Streamwise maximum detection

### Stored-grid maximum

```matlab
cfg.maxDetection.method = 'grid';
```

This checks only stored streamwise stations and generally provides a lower approximation to continuous first-transition probability.

### Cubic-Hermite maximum

```matlab
cfg.maxDetection.method = 'hermite';
```

This uses endpoint values and derivatives of each directional gain to locate interior maxima between stations.

For precomputed OWNS data, derivatives are currently estimated with a nonuniform three-point polynomial stencil.

### Adaptive method

```matlab
cfg.maxDetection.method = 'adaptive';
```

This branch is reserved for future work and is not currently part of the verified production workflow.

## Angular integration

### Two-dimensional deterministic rule

For \(r=2\):

```matlab
cfg.angular.method = 'circle';
cfg.angular.numDirections = 4096;
```

### Scrambled Sobol’ RQMC

For moderate or high stochastic dimension:

```matlab
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 4096;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 1;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;
```

RQMC replicates are evaluated using common local and first-transition directional samples, preserving:

\[
p_{\mathrm{local}}(x)
\leq
F_{X_{\mathrm{tr}}}(x)
\]

replicate by replicate.

## Transition-oriented dimension reduction

The full OWNS stochastic rank is currently \(r=141\). A fixed reduced basis is used over the complete streamwise domain.

### Integrated basis

```matlab
cfg.reduction.method = 'integrated';
cfg.reduction.targetRank = 30;
cfg.reduction.integrationWeight = 'uniform';
```

The integrated matrix is:

\[
K
=
\int G(x)\,dx.
\]

Its leading eigenvectors define the reduced stochastic basis.

### Other available basis methods

```text
eigUnion
hybrid
greedy
givenBasis
none
```

For the current real-coefficient OWNS dataset, a coupled RQMC rank study indicated that rank 30 satisfies an approximate absolute complete-CDF tolerance of \(5\times10^{-3}\) under the provisional threshold.

Rank adequacy should be checked again if:

- the dataset changes materially;
- the threshold becomes streamwise varying;
- the energy norm changes;
- or the Gaussian coefficient convention changes.

## Local energy probabilities

Two local probability methods are available.

### Same-angular-rule estimate

This uses the same RQMC directions as the first-transition calculation:

```matlab
cfg.localProbability.method = 'sameAngular';
```

### Independent characteristic-function inversion

```matlab
cfg.localProbability.method = 'characteristicFunction';
```

or:

```matlab
cfg.localProbability.method = 'both';
```

The deterministic method evaluates the local generalized quadratic-form tail using characteristic-function inversion and is independent of angular RQMC.

## Threshold and amplitude studies

For a constant threshold and inlet amplitude multiplier \(\epsilon\),

\[
F_{X_{\mathrm{tr}}}(x)
=
\mathbb E_u
\left[
\overline F_{\chi_r^2}
\left(
\frac{e_{\mathrm{thres}}}
{\epsilon^2 M_u(x)}
\right)
\right].
\]

A reusable unnormalized energy envelope permits many thresholds and amplitudes to be evaluated without repeating the directional maximum calculation.

Enable the main-script sweep with:

```matlab
cfg.thresholdSweep.enable = true;
```

Representative options are:

```matlab
cfg.thresholdSweep.factors = logspace(-0.5, 0.5, 9);
cfg.thresholdSweep.amplitudes = [0.5, 0.75, 1.0, 1.25, 1.5];
cfg.thresholdSweep.referenceMethod = 'maximumMeanEnergy';
```

A target terminal probability can also be used to solve for a constant threshold.

## Tests

Run the following tests from the MATLAB command window.

### Synthetic regression

```matlab
runSyntheticRegressionTests
```

Verifies:

- construction of the synthetic \(G(x)\);
- angular convergence;
- trajectory Monte Carlo agreement;
- local-versus-first-transition consistency;
- probability conservation.

### Hermite maximum detection

```matlab
runHermiteMaximumTests
```

Verifies:

- exact reproduction of cubic maxima;
- exact synthetic derivatives;
- Hermite versus stored-grid maxima;
- agreement with a dense streamwise reference;
- finite-difference derivative fallback.

### RQMC integration

```matlab
runRQMCTests
```

Verifies:

- scrambled Sobol’ convergence;
- deterministic reproducibility under fixed seeds;
- agreement with a deterministic circular reference;
- replicate uncertainty.

### OWNS smoke test

```matlab
runOWNSSmokeTest
```

Verifies the complete OWNS loading, energy compression, threshold, Hermite, and RQMC path.

### Real-model rank convergence

```matlab
runOWNSCoupledReductionStudy
```

Uses common RQMC coordinates across nested ranks to estimate paired reduction error.

### Proper-complex verification

```matlab
runProperComplexRegressionTests
```

Verifies:

- analytical scalar proper-complex probabilities;
- proper-complex energy moments;
- direct trajectory Monte Carlo agreement.

### Proper-complex rank convergence

```matlab
runOWNSProperComplexRankStudy
```

Checks complex-specific transition-oriented bases and full-rank convergence.

### Generalized quadratic-form tails

```matlab
runQuadraticFormTailTests
```

Verifies deterministic local tail probabilities and threshold inversion.

### OWNS local probability comparison

```matlab
runOWNSLocalProbabilityCheck
```

Compares angular RQMC local probabilities against independent characteristic-function inversion.

### Threshold sensitivity

```matlab
runOWNSThresholdStudy
```

Evaluates constant-threshold and inlet-amplitude sensitivity using a reusable energy envelope.

## Representative verified results

### Synthetic two-dimensional test

Endpoint-only circular result:

```text
F(xmax) = 0.3461435512
```

Hermite-corrected result:

```text
F(xmax) = 0.3461875
```

### Current OWNS dataset with provisional threshold

Using the real-coefficient model and the provisional automatic threshold, the full-rank transition probability was approximately:

```text
F(xmax) ≈ 0.447
```

Using proper-complex coefficients produced a larger probability of approximately:

```text
F(xmax) ≈ 0.483
```

These values are for numerical development only because the physical transition threshold has not yet been calibrated.

## Current limitations

The following capabilities are not currently part of the verified production path:

- direct execution of the OWNS spatial march;
- nonzero distributed forcing;
- adaptive insertion of new physical OWNS stations;
- nonlinear disturbance propagation;
- DNS threshold calibration;
- broadband joint reconstruction across multiple \((\omega,\beta)\) pairs;
- directory-wide processing of all spectral pairs;
- physical-field plotting and reconstruction inside this repository;
- automatic parallel RQMC execution.

OWNS runs and physical disturbance visualization are expected to remain external to this codebase.

## Planned future work

Current priorities are:

1. finish verification of deterministic generalized quadratic-form tails;
2. add configurable higher-order derivative stencils for \(G'(x)\);
3. automate RQMC convergence for production calculations;
4. add optional parallel execution over RQMC replicates;
5. implement survivor-conditioned stochastic covariance;
6. implement transition-conditioned directional diagnostics;
7. design directory-wide processing of all \((\omega,\beta)\) pairs;
8. design joint reconstruction over the full spectral space.

The broadband reconstruction requires a deliberate modeling decision about spectral correlations, phase conventions, and whether transition is defined from summed energy or another joint spectral statistic.

## Reproducibility

For reproducible calculations, preserve:

- the complete configuration structure;
- selected OWNS file path;
- temporal frequency and spanwise wavenumber;
- stochastic rank and basis;
- threshold definition;
- maximum-detection method;
- RQMC point count;
- RQMC replicate count;
- all random seeds;
- MATLAB version;
- numerical tolerances.

The main result structure should include, at minimum:

```matlab
results.cfg
results.meta
results.prob
results.post
results.reductionInfo
```

## Data and repository policy

Large OWNS MAT-files should not be committed to Git. Store them externally and configure their paths at runtime.

Likewise, generated result files, figures, and temporary MATLAB files should generally be excluded unless intentionally added as small regression references.
```

```gitignore
# MATLAB temporary and autosave files
*.asv
*.m~
*.autosave

# MATLAB generated files
*.mex*
*.mlappinstall
*.mltbx

# Large data and result files
*.mat
*.h5
*.hdf5

# Generated output directories
results/
figures/
output/
tmp/
temp/

# Keep empty project data directory
!data/.gitkeep

# Profiling and coverage output
profile_results/
coverage/

# Operating-system files
.DS_Store
Thumbs.db

# Editor and IDE settings
.vscode/
.idea/

# Slurm or cluster output
slurm-*.out
*.o
*.e
```

If you intend to commit small MAT-files for regression tests, replace the broad `*.mat` rule with:

```gitignore
data/*.mat
results/*.mat
output/*.mat

!tests/reference_data/*.mat
