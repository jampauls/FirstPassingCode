% OWNS_run_first_transition_ensemble.m
%
% Driver for the total-energy first-transition probability across an
% ensemble of (omega,beta) OWNS solution files.
%
% Each mode file is treated as an independent stochastic contribution:
%
%   e_total(x) = sum_i w_i' A_i(x) w_i,     w_i ~ N(0,I_{r_i}) independent
%
% The first-transition probability is defined from e_total(x) crossing a
% single threshold, rather than from any individual mode's energy.
%
% MATLAB version target: R2020b

clear; %clc;

%% ========================================================================
%  0. Path setup
% ========================================================================

addpath(genpath(fullfile(pwd, 'src')));
addpath(genpath(fullfile(pwd, 'examples')));
addpath(genpath(fullfile(pwd, 'data')));
addpath(genpath(fullfile(pwd, 'tests')));
addpath(genpath(fullfile(pwd, 'figures')));

%% ========================================================================
%  1. User configuration
% ========================================================================

cfg = struct();

% -------------------------------------------------------------------------
% Ensemble file discovery
% -------------------------------------------------------------------------

cfg.owns = struct();

% All MAT-files in this directory are treated as ensemble mode files,
% regardless of filename.
cfg.owns.ensembleDirectory = fullfile('data');

% Mode whose streamwise grid defines the common coordinate (index into
% the sorted file listing, or an explicit file path). Defaults to the
% first file if left empty.
cfg.owns.ensembleReferenceFile = [];

% Directory holding the cached minimal {x, A, meta} form of each raw mode
% file. Built automatically on first access to each raw file and reused
% on every later run, avoiding repeated multi-GB reads of the raw OWNS
% solutions.
cfg.owns.reducedCacheDirectory = fullfile('data', 'reducedCache');

% Target number of streamwise stations retained when building the cache.
% The largest OWNS files carry ~10,000 stations to resolve the march
% itself, far more than needed to sample the propagated stochastic
% structure for first-transition prediction. Leave empty to keep every
% station (no downsampling).
cfg.owns.numCacheStations = 400;

% Optional 1/sqrt(omega) domain cutoff applied before the station subset
% above is chosen: only stations with (x - x(1)) <= coefficient/sqrt(|w|)
% are eligible. Useful for concentrating a fixed station budget near the
% inlet for high-frequency modes, whose relevant dynamics occur over a
% correspondingly shorter streamwise distance. Leave empty to use the
% full streamwise domain regardless of omega.
cfg.owns.cacheXCutoffCoefficient = [];

cfg.owns.solutionVariable = 'solution';
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;
cfg.owns.validateWeights = true;
cfg.owns.weightTolerance = 1e-12;
cfg.owns.verbose = false;

% -------------------------------------------------------------------------
% Streamwise-coordinate configuration
% -------------------------------------------------------------------------

cfg.coordinates = struct();
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.referenceIndex = 1;
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;

% -------------------------------------------------------------------------
% Total-energy transition-threshold configuration
% -------------------------------------------------------------------------

cfg.threshold = struct();

% Options:
%   'specifiedScalar'     : cfg.threshold.value, an absolute threshold
%   'initialMeanMultiple' : cfg.threshold.multiplier * inlet mean total energy
%   'auto'                : cfg.threshold.autoMaxMeanFactor * max mean total energy
cfg.threshold.method = 'auto';
cfg.threshold.value = [];
cfg.threshold.multiplier = 10;
cfg.threshold.autoMaxMeanFactor = 1;
cfg.threshold.warnIfProvisional = true;

% -------------------------------------------------------------------------
% Per-mode dimension reduction
% -------------------------------------------------------------------------

cfg.reduction = struct();

% Reduction method applied independently to every mode.
cfg.reduction.method = 'integrated';

% Reduced rank kept per mode. The total stochastic dimension is
% numModes * perModeTargetRank, and must not exceed the Sobol'
% construction limit (1111) when cfg.angular.method = 'rqmc'.
cfg.reduction.perModeTargetRank = 4;

cfg.reduction.integrationWeight = 'uniform';
cfg.reduction.numIntegratedModes = 10;
cfg.reduction.numModesPerStation = 3;
cfg.reduction.stationSelection = 'recordsAndUniform';
cfg.reduction.numUniformStations = 12;
cfg.reduction.recordRelativeTolerance = 1e-6;
cfg.reduction.useGreedyEnrichment = true;
cfg.reduction.greedyTolerance = 1e-8;
cfg.reduction.maxGreedyModes = 50;
cfg.reduction.snapshotTolerance = 1e-12;
cfg.reduction.saveBasis = false;
cfg.reduction.basisFile = '';
cfg.reduction.verbose = false;

% -------------------------------------------------------------------------
% Streamwise maximum detection
% -------------------------------------------------------------------------

cfg.maxDetection.method = 'hermite';
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;

% -------------------------------------------------------------------------
% Angular integration
% -------------------------------------------------------------------------

cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 4096;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 1;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;
cfg.angular.autoConverge = false;

% -------------------------------------------------------------------------
% Numerical options
% -------------------------------------------------------------------------

cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;
cfg.numerics.useUpperTail = true;

% -------------------------------------------------------------------------
% Output and plotting
% -------------------------------------------------------------------------

cfg.plot.enable = true;
cfg.plot.showLocalProbability = true;
cfg.plot.showIntervalProbabilities = true;

%% ========================================================================
%  2. Run the ensemble workflow
% ========================================================================

results = runOWNSEnsembleWorkflow(cfg);

x = results.x;
prob = results.prob;

fprintf('\nEnsemble summary:\n');
fprintf('  modes combined        = %d\n', results.numModes);
fprintf('  total stochastic dim  = %d (full rank %d)\n', ...
    results.rTotal, results.rFull);
fprintf('  threshold method      = %s\n', results.thresholdInfo.method);
fprintf('  total threshold value = %.6e\n', results.thresholdInfo.maximumThreshold);

%% ========================================================================
%  3. Post-processing and plotting
% ========================================================================

post = struct();
post.quantileLevels = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
post.quantiles = computeTransitionQuantiles(x, prob.F, post.quantileLevels);

pLocal = [];
post.memoryCorrection = [];

if isfield(prob, 'pLocal')
    pLocal = prob.pLocal;
    post.memoryCorrection = prob.F(:) - pLocal(:);
end

if cfg.plot.enable
    plotTransitionResults(x, prob, pLocal, post, cfg);
end
