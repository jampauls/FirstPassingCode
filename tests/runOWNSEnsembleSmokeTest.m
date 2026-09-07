% runOWNSEnsembleSmokeTest.m
%
% End-to-end test of the multi-mode total-energy first-transition
% workflow (runOWNSEnsembleWorkflow), using every OWNS solution file
% currently in data/ as an independent (omega,beta) mode.
%
% Also exercises the reduced-cache layer in loadOWNSModeCached: the
% workflow is run twice against a fresh cache directory, and the second
% run is checked to reuse the cached files rather than rebuilding them.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS ensemble first-transition smoke test\n');
fprintf('============================================================\n\n');

dataDir = fullfile(projectRoot, 'data');
cacheDir = fullfile(tempdir, 'FirstPassingCodeEnsembleSmokeTestCache');

if isfolder(cacheDir)
    rmdir(cacheDir, 's');
end

%% Configuration

cfg = struct();

cfg.owns = struct();
cfg.owns.ensembleDirectory = dataDir;
cfg.owns.ensembleReferenceFile = [];
cfg.owns.reducedCacheDirectory = cacheDir;
cfg.owns.solutionVariable = 'solution';
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;
cfg.owns.validateWeights = true;
cfg.owns.weightTolerance = 1e-12;
cfg.owns.verbose = false;

cfg.coordinates = struct();
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.referenceIndex = 1;
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;

cfg.threshold = struct();
cfg.threshold.method = 'auto';
cfg.threshold.value = [];
cfg.threshold.multiplier = 10;
cfg.threshold.autoMaxMeanFactor = 1;
cfg.threshold.warnIfProvisional = true;

cfg.reduction = struct();
cfg.reduction.method = 'integrated';
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

cfg.maxDetection = struct();
cfg.maxDetection.method = 'hermite';
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;

cfg.angular = struct();
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 256;
cfg.angular.numReplicates = 4;
cfg.angular.randomSeed = 1234;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;
cfg.angular.autoConverge = false;

cfg.numerics = struct();
cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;
cfg.numerics.useUpperTail = true;

%% Expected input files

fileList = discoverOWNSEnsembleFiles(dataDir);

fprintf('Found %d ensemble mode files in %s:\n', numel(fileList), dataDir);
for k = 1:numel(fileList)
    [~, fname, fext] = fileparts(fileList{k});
    fprintf('  %s%s\n', fname, fext);
end
fprintf('\n');

if numel(fileList) < 2
    error('runOWNSEnsembleSmokeTest:TooFewFiles', ...
        'Expected at least 2 OWNS mode files in data/, found %d.', ...
        numel(fileList));
end

%% First run: builds the reduced cache from raw files

fprintf('--- First run (cold cache) ---\n\n');

resultsFirst = runOWNSEnsembleWorkflow(cfg);

validateProbabilityCurves(resultsFirst.prob, cfg);

if resultsFirst.numModes ~= numel(fileList)
    error('runOWNSEnsembleSmokeTest:ModeCountMismatch', ...
        'Expected %d modes, workflow combined %d.', ...
        numel(fileList), resultsFirst.numModes);
end

expectedRTotal = resultsFirst.numModes * cfg.reduction.perModeTargetRank;
if resultsFirst.rTotal > expectedRTotal
    error('runOWNSEnsembleSmokeTest:UnexpectedRank', ...
        'Total stochastic dimension %d exceeds %d modes * rank %d.', ...
        resultsFirst.rTotal, resultsFirst.numModes, ...
        cfg.reduction.perModeTargetRank);
end

cacheListing = dir(fullfile(cacheDir, '*.mat'));
if numel(cacheListing) ~= numel(fileList)
    error('runOWNSEnsembleSmokeTest:CacheNotBuilt', ...
        'Expected %d cached mode files, found %d in %s.', ...
        numel(fileList), numel(cacheListing), cacheDir);
end

fprintf('Reduced cache populated with %d files in:\n  %s\n\n', ...
    numel(cacheListing), cacheDir);

%% Second run: must reuse the cache, not rebuild it

fprintf('--- Second run (warm cache) ---\n\n');

cacheInfoBefore = arrayfun(@(f) f.datenum, cacheListing);

resultsSecond = runOWNSEnsembleWorkflow(cfg);

cacheListingAfter = dir(fullfile(cacheDir, '*.mat'));
cacheInfoAfter = arrayfun(@(f) f.datenum, cacheListingAfter);

if numel(cacheListingAfter) ~= numel(cacheListing) || ...
        any(sort(cacheInfoAfter) ~= sort(cacheInfoBefore))
    error('runOWNSEnsembleSmokeTest:CacheRebuilt', ...
        'Cached mode files were rewritten on the second run.');
end

validateProbabilityCurves(resultsSecond.prob, cfg);

%% Cross-check: cached run must reproduce the cold-cache result

if ~isequal(size(resultsFirst.prob.F), size(resultsSecond.prob.F))
    error('runOWNSEnsembleSmokeTest:SizeMismatch', ...
        'First- and second-run probability curves have different sizes.');
end

cdfDefect = max(abs(resultsFirst.prob.F - resultsSecond.prob.F));

if cdfDefect > 1e-12
    error('runOWNSEnsembleSmokeTest:CacheChangedResult', ...
        'Cached run differs from cold-cache run by %.3e.', cdfDefect);
end

%% Report and clean up

fprintf('Ensemble smoke-test result:\n');
fprintf('  modes combined        = %d\n', resultsFirst.numModes);
fprintf('  total stochastic dim  = %d (full rank %d)\n', ...
    resultsFirst.rTotal, resultsFirst.rFull);
fprintf('  threshold method      = %s\n', resultsFirst.thresholdInfo.method);
fprintf('  F(xmax)               = %.12e\n', resultsFirst.prob.F(end));
fprintf('  S(xmax)               = %.12e\n', resultsFirst.prob.S(end));
fprintf('  cache reuse defect    = %.3e\n\n', cdfDefect);

rmdir(cacheDir, 's');

fprintf('PASSED: runOWNSEnsembleSmokeTest\n');
