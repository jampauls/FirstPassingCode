# OWNS First-Transition Prediction

MATLAB R2020b workflow for estimating the streamwise location of the first energy-threshold crossing in a propagated stochastic OWNS solution. The main entry point is [exact_Struct_run_first_transition_latest.m](exact_Struct_run_first_transition_latest.m).

For a real Gaussian stochastic factor $w \sim N(0, I_r)$, the script forms

$$
q(x) = B(x)w, \qquad e(x) = q(x)^* H(x)q(x), \qquad
G(x) = \frac{\operatorname{real}(B(x)^*H(x)B(x))}{e_{\mathrm{thres}}(x)}.
$$

It evaluates the memory-aware first-transition CDF from directional running maxima:

$$
F_{X_{tr}}(x) = \mathbb{E}_u\left[\chi_r^2\!\operatorname{sf}\left(
\frac{1}{\max_{\xi \leq x} u^T G(\xi)u}\right)\right].
$$

## Requirements

- MATLAB R2020b or compatible release.
- Statistics and Machine Learning Toolbox for `sobolset` and chi-square tail evaluation when using the default RQMC calculation.
- An OWNS solution MAT-file containing a variable named `solution`.

The checked-in default input is [data/ninth_Run_i_1_j_1_solution.mat](data/ninth_Run_i_1_j_1_solution.mat). To use another data set, change `cfg.dataFile` near the top of the main script. The script uses `cfg.owns.dataFile`, initialized from that setting, for the default `ownsmatfile` source.

## Input Data Contract

The `solution` structure must contain:

| Field | Expected role |
| --- | --- |
| `q` | Propagated factor with size `[6*Ny, Nx, r]`. |
| `q_in` | Inlet energy-bearing factor with size `[5*Ny, r]`. |
| `gram_W` | Nonnegative energy weights with size `[5*Ny, Nx]`. |
| `x`, `y` | Physical grid arrays with size `[Ny, Nx]`. |
| `eta` | Wall-normal grid vector of length `Ny`. |
| `w` | Temporal frequency. |
| `delta` | Required when delta99 trimming is enabled; one positive delta99 value per streamwise station. |
| `beta` | Optional spanwise wavenumber. |
| `z` | Optional grid used by arc-length coordinates. |

The energy uses the first five state blocks `[rho, u, v, w, T]`; the pressure block is excluded. The loader validates dimensions and energy-weight signs, converts the solution into compressed stochastic energy matrices $A(x)$, and does not retain the large physical factor after loading.

## Run

From the repository root in MATLAB:

```matlab
run('exact_Struct_run_first_transition_latest.m')
```

The current default configuration:

- loads the OWNS MAT-file (`cfg.dataSource = 'ownsmatfile'`),
- uses wall arc length as the streamwise coordinate,
- creates an automatic constant threshold equal to the maximum mean energy,
- reduces the stochastic dimension with a rank-30 integrated basis,
- finds continuous interval maxima using cubic Hermite interpolation, and
- estimates the angular integral with 16 scrambled Sobol replicates of 8192 directions each.

The default calculation can be expensive. For a smaller end-to-end check, run:

```matlab
run('tests/runOWNSSmokeTest.m')
```

That test uses 4 RQMC replicates of 512 directions, validates the probability curves, and plots basic diagnostics. It has its own `cfg.owns.dataFile`, which may need to be updated for a local OWNS file.

## Configuration

Edit the `cfg` structure in the main script before running.

| Setting | Main choices |
| --- | --- |
| `cfg.dataSource` | `ownsmatfile` (default), `ownsUserSelect`, `synthetic`, or `matfile`. `ownsmarch` is not implemented. |
| `cfg.maxDetection.method` | `grid`, `hermite` (default), or `adaptive`. |
| `cfg.angular.method` | `rqmc` (default), `random`, `circle` (for `r = 2`), or `user`. |
| `cfg.reduction.method` | `none`, `givenBasis`, `integrated` (default), `eigUnion`, `hybrid`, or `greedy`. |
| `cfg.threshold.method` | `auto` (default), scalar/vector-specified, mean-energy-based, or inlet-quantile approximation. |
| `cfg.owns.trimToDelta99` | Enables an energy-domain restriction to a near-wall region; default `false`. |
| `cfg.owns.numDelta99` | Positive number of local delta99 thicknesses retained when trimming is enabled; default `3`. |
| `cfg.owns.delta99MaskMethod` | `cosineTaper` (default in the driver) or `binary` for an exact cutoff. |
| `cfg.owns.delta99TaperWidth` | Half-width of the cosine transition band in delta99 units; default `0.25`. |

With `cosineTaper`, weights are one below
$\texttt{numDelta99}\,\delta_{99}(x)-\texttt{delta99TaperWidth}\,\delta_{99}(x)$,
decay with a raised cosine through the transition band, and are zero above
the outer edge. `binary` sets weights outside
$|y-y_{wall}| \leq \texttt{numDelta99}\,\delta_{99}(x)$ to zero exactly.
Both options restrict the energy inner product without changing the
propagated state or stochastic coordinates.

Set `cfg.output.saveResults = true` to save a `results` structure to `cfg.output.resultsFile` (default: `data/results_first_transition.mat`). It includes the effective configuration, coordinates, probability curves, post-processing, source metadata, angular metadata, maximum-method metadata, reduction information, and the reduction basis when one is used.

## Results

The driver reports:

- terminal first-transition probability `prob.F(end)`,
- right-censoring probability `prob.S(end)`,
- terminal RQMC standard error for replicate-mean RQMC,
- local-exceedance consistency against the first-transition CDF, and
- transition quantiles and directional diagnostics in `post`.

With plotting enabled, it calls `plotTransitionResults`. For synthetic data, it also calls `plotSyntheticProblem`.

## Workflow Components

The main script depends on the following project areas:

| Directory | Role in this workflow |
| --- | --- |
| `src/data` | OWNS loading, coordinate extraction, threshold construction, and input validation. |
| `src/matrices` | Compressed energy matrices, normalized matrix family construction, and derivatives. |
| `src/reduction` | Transition-oriented stochastic basis construction and matrix/factor projection. |
| `src/angular` | Random, deterministic, and scrambled-Sobol sphere directions. |
| `src/maxima` | Grid, Hermite, and adaptive directional running maxima. |
| `src/probability` | RQMC CDF estimation, survival/CDF construction, local exceedance, and quantiles. |
| `src/diagnostics` | Directional post-processing diagnostics. |
| `src/plotting` | Transition-result visualizations. |

Other scripts and studies in the repository are intentionally outside this README's scope.
