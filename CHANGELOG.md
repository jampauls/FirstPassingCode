# Changelog

## Ensemble branch (this session)

### Added

- **Multi-mode ensemble first-transition workflow** (`OWNS_run_first_transition_ensemble.m`, `src/ensemble/runOWNSEnsembleWorkflow.m`): combines an arbitrary number of independent \((\omega,\beta)\) OWNS solution files and computes the first-transition probability for the *total* disturbance energy summed across all modes, instead of any single mode's energy. Modes are combined as a block-diagonal quadratic form and reused unchanged through the existing angular/maxima/probability machinery.
  - `src/ensemble/discoverOWNSEnsembleFiles.m` — lists all mode MAT-files in a directory.
  - `src/ensemble/loadOWNSModeRaw.m` — loads one mode and compresses it to `{x, A, meta}` without a per-mode threshold.
  - `src/ensemble/buildEnsembleTransitionThreshold.m` — defines a threshold from the total (summed) mean/variance energy.
  - `src/ensemble/interpMatrixFamilyToGrid.m` — interpolates each mode's energy matrices onto a common streamwise grid.
  - `src/ensemble/reduceEnsembleMode.m` — applies per-mode transition-oriented stochastic dimension reduction.
  - `src/ensemble/assembleBlockDiagonalFamily.m` — assembles reduced per-mode matrices into a block-diagonal family.
  - Fixed `src/matrices/buildGprimeFamily.m` to accumulate with a sparse zero matrix when `G` is sparse, preserving block-diagonal sparsity for the ensemble path.

- **Persistent, staleness-aware reduced-data cache** (`src/ensemble/loadOWNSModeCached.m`): builds the small `{x, A, meta}` form of each raw mode file once and caches it under `cfg.owns.reducedCacheDirectory` (default `data/reducedCache`), avoiding repeated multi-gigabyte reads of raw OWNS files across passes and across runs. The cache is automatically invalidated and rebuilt if a raw file's modification time or size changes; new files are cached without disturbing existing entries.

- **Tunable streamwise-station downsampling** (`src/ensemble/downsampleOWNSStations.m`): reduces the retained station count before any per-station energy matrix is built and before the cache is saved, via `cfg.owns.numCacheStations`. An optional `cfg.owns.cacheXCutoffCoefficient` restricts eligible stations to `(x - x(1)) <= coefficient / sqrt(|omega|)`, concentrating resolution near the inlet for high-frequency modes.

- **Tests**:
  - `tests/runOWNSEnsembleSmokeTest.m` — end-to-end ensemble workflow test against every file in `data/`; verifies cold/warm cache reuse and identical results.
  - `tests/runOWNSCacheStalenessTest.m` — verifies cache persistence, staleness detection/rebuild on source-file changes, and isolated caching of newly added files.
  - `tests/runOWNSCacheDownsamplingSweepTest.m` — sweeps 5 station-count resolutions, verifying full-domain span at low resolution and omega-scaled domain truncation at high resolution.

### Fixed

- `downsampleOWNSStations.m`: the "no downsampling needed" shortcut previously only compared `numCacheStations` to the raw station count, silently ignoring the `cacheXCutoffCoefficient` cutoff whenever the requested station count exceeded the raw grid size.

### Documentation

- Updated `README.md` with the multi-mode ensemble workflow, reduced-data cache, and station-downsampling sections, and refreshed the repository structure, test list, capabilities, limitations, and planned-work sections accordingly.
