# Research Memo: OWNS First-Transition Probability Codebase

**Purpose.** This memo describes the mathematical method, software design, and performance features of the OWNS first-transition probability codebase. It excludes problem-specific configuration and numerical results.

## Executive Summary

The code estimates the probability that disturbance energy first crosses a specified threshold at, or before, each streamwise location. It uses a low-dimensional Gaussian representation of propagated OWNS disturbances. The method preserves upstream crossing history. Thus, it differs from a local exceedance calculation at each station.

The main strategy separates Gaussian magnitude from Gaussian direction. The radial probability is evaluated analytically with a chi-square tail. Numerical integration remains only over unit directions. The code supports deterministic, random, and randomized quasi-Monte Carlo (RQMC) angular rules.

The code reduces cost through stochastic dimension reduction, streamwise station downsampling, persistent reduced caches, cubic-Hermite maximum resolution, and bounded parallel cache construction. The primary multi-mode driver is `OWNS_run_first_transition_ensemble.m`.

## Mathematical Formulation

### Gaussian disturbance model

At streamwise location $x$, the propagated disturbance is represented as

$$
q(x) = B(x)w, \qquad w \sim N(0, I_r).
$$

The default convention uses real Gaussian coordinates $w$, even when OWNS factors and state vectors are complex. Given energy weight matrix $H(x)$, the disturbance energy becomes

$$
e(x) = q(x)^* H(x) q(x) = w^T A(x)w,
$$

where

$$
A(x) = \operatorname{Re}\{B(x)^*H(x)B(x)\}.
$$

For threshold $e_{\mathrm{thres}}(x)>0$, the normalized energy matrix is

$$
G(x) = \frac{A(x)}{e_{\mathrm{thres}}(x)}.
$$

The first-transition location is

$$
X_{\mathrm{tr}} = \inf\{x : w^T G(x)w \geq 1\}.
$$

### Directional reduction

The Gaussian vector is decomposed into radius and direction:

$$
w = Ru, \qquad R^2 \sim \chi_r^2, \qquad u \sim \operatorname{Unif}(S^{r-1}).
$$

For each direction, define the history-aware directional maximum

$$
m_u(x) = \max_{\xi \leq x} u^T G(\xi)u.
$$

The first-transition CDF is

$$
F_{X_{\mathrm{tr}}}(x) = E_u\left[\overline F_{\chi_r^2}\left(\frac{1}{m_u(x)}\right)\right].
$$

This is the central result used by the code. It converts a trajectory event in Gaussian space into an angular expectation. A direction that has crossed the threshold remains transitioned downstream because $m_u(x)$ is nondecreasing.

The code also computes survival probability $S(x)=1-F(x)$, interval transition probability, right-censoring probability, discrete hazard, and local exceedance probability. The local value omits the running maximum and must not exceed the first-transition CDF when both use the same samples.

### Multi-mode ensemble model

For independent OWNS modes, total energy is the sum of mode energies:

$$
e_{\mathrm{total}}(x) = \sum_{i=1}^{M} w_i^T A_i(x)w_i.
$$

The code reduces each mode independently, then forms a block-diagonal normalized matrix family. This preserves mode independence while allowing the existing maximum and probability routines to operate on a combined stochastic vector.

## Conditional Monte Carlo and Novelty

Direct trajectory Monte Carlo estimates $F_{\mathrm{fp}}(x)$ with the indicator

$$
Y_x=\mathbf{1}\left\{
\sup_{\xi\leq x}\boldsymbol{w}^TG(\xi)\boldsymbol{w}\geq1
\right\}.
$$

Conditioning on the Gaussian direction $\boldsymbol{u}$ analytically integrates the radius:

$$
E[Y_x\mid\boldsymbol{u}]=
\overline F_{\chi_r^2}\left(1/m_{\boldsymbol{u}}(x)\right).
$$

Thus the angular integrand is the conditional Monte Carlo, or Rao-Blackwellized, form of the direct indicator. The tower property preserves the target probability. The law of total variance gives variance no greater than direct independent trajectory Monte Carlo. This use of Rao-Blackwellization does not require $\boldsymbol{u}$ to be a sufficient statistic for an unknown parameter.

The established part is the Gaussian radial-angular decomposition and conditional-expectation variance reduction. The method-specific contribution is their explicit combination with the running maximum of a spatial family of quadratic forms. This construction preserves first-passage history while removing the radial random variable analytically. The same form extends to independent multi-mode systems through block-diagonal quadratic matrices, and to proper-complex Gaussian coordinates through a gamma radial tail.

## Implementation Structure

The workflow has six main stages.

1. Discover mode files and prepare a reduced cache for each source file.
2. Select one mode as the common streamwise reference grid.
3. Interpolate mode energy matrices to that grid and accumulate mean and variance energy for threshold construction.
4. Normalize each mode by the common threshold, reduce its stochastic rank, and assemble the block-diagonal matrix family.
5. Build $G'(x)$ from supplied derivatives or three-point finite differences when Hermite or adaptive maxima are requested.
6. Evaluate directional running maxima and integrate chi-square tail probabilities over directions.

The workflow returns probability curves, reduction metadata, energy diagnostics, and stage timings in `results.timings`. It validates matrix dimensions, angular weights, probability normalization, monotonicity, and the condition that a local gain does not exceed its running maximum.

## Numerical Methods

### Angular integration

For two stochastic dimensions, the code can use an equally spaced projective-circle rule. Higher-dimensional calculations can use independent Gaussian-normalized directions or scrambled Sobol' RQMC directions. RQMC uses independently scrambled replicates. Their mean is the estimate; their sample variation gives a standard-error estimate.

Each RQMC replicate has a deterministic seed derived from the base seed and replicate index. This makes repeated serial runs reproducible. The Sobol' construction limits RQMC to 1111 stochastic dimensions.

### Streamwise maxima

The grid method evaluates maxima only at saved stations. The Hermite method uses endpoint gains and directional derivatives to construct a cubic on each streamwise interval. It evaluates interval endpoints and internal stationary points. This resolves an interpolated maximum without globally refining the streamwise grid.

The running maximum update is sequential in streamwise position. It is not a direct `parfor` target because each station depends on the previous maximum.

### Stochastic dimension reduction

The reduction subsystem selects a basis that retains transition-relevant energy structure. Supported strategies include integrated, eigenspace-union, hybrid, and greedy bases. Projection reduces matrix size before angular integration. In ensemble work, this occurs independently for each mode.

## Performance Methods

### Persistent reduced cache

Reading a raw OWNS file and constructing $A(x)$ can dominate runtime and memory use. The cache stores only the streamwise coordinate, energy matrices, and required metadata. It is invalidated when the source file timestamp or size changes.

Cache publication is hardened against interrupted writes. The code saves to a temporary file in the cache directory, then moves it to the final name. An unreadable or incomplete cache file is ignored and rebuilt. This prevents a partially written final cache file from being consumed by later stages.

### Parallel cache preparation

Before the serial numerical workflow, independent missing or stale cache entries can be built with `parfor`. The worker count is bounded by the number of files, available CPU cores, available memory, and an optional user cap. It reserves two CPU cores and at least 25% of physical memory or 4 GB. The memory estimate allocates three times the largest raw file size per worker. If one worker is the safe limit, or parallel execution is unavailable, the same operation runs serially.

This parallel phase is appropriate because each mode has a distinct cache file. Other dominant stages remain serial because they aggregate data or have streamwise dependencies.

### Station downsampling and domain restriction

The cache builder can retain a limited number of streamwise stations before it forms energy matrices. This reduces both matrix-construction work and cache size. An optional cutoff proportional to $1/\sqrt{|\omega|}$ restricts the eligible domain before station selection. The first and final eligible stations are retained.

### Batching and stable probability evaluation

Directional gains use matrix-matrix products, which batch many directions at a station. The probability routines evaluate both chi-square tails with `gammainc` and accumulate the smaller tail directly. This reduces loss of precision when a probability is close to zero or one.

### Timing instrumentation

The ensemble workflow measures cache preparation, reference setup, pass 1, reduction and assembly, derivative construction, and angular integration. The timings support focused profiling without changing the numerical method.

## Limits and Development Priorities

The reduced cache trades raw OWNS resolution for a selected streamwise representation. The station budget and reduction rank require convergence checks for a new problem class. Hermite maxima improve interval resolution, but use interpolation rather than a certified bound. RQMC uncertainty is estimated from replicates and remains subject to finite-sample error.

Timing output should guide further work. In a cold-cache run, cache creation is often the main target. In a warm-cache run, reduction and angular integration can become more important. Any future parallelization should keep deterministic RQMC seeds and should not break the sequential running maximum definition.