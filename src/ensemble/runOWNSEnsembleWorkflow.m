function results = runOWNSEnsembleWorkflow(cfg)
%RUNOWNSENSEMBLEWORKFLOW First-transition probability for many OWNS modes.
%
%   results = runOWNSEnsembleWorkflow(cfg)
%
% Combines many (omega,beta) OWNS solution files, each treated as an
% independent stochastic mode, into the total disturbance energy
%
%   e_total(x) = sum_i w_i' A_i(x) w_i,     w_i ~ N(0,I_{r_i}) independent.
%
% Because the modes are independent, e_total(x) is a quadratic form in
% the concatenated vector [w_1;...;w_M] with block-diagonal matrix
% A_total(x) = blkdiag(A_1(x),...,A_M(x)). The existing angular/maxima/
% probability machinery is reused unchanged by keeping this block
% structure sparse.
%
% Required cfg fields:
%   cfg.owns.ensembleDirectory   directory containing one MAT-file per mode
%   cfg.owns.ensembleReferenceFile (optional) index or path selecting the
%                                  mode whose grid defines the common
%                                  streamwise coordinate
%   cfg.coordinates               as in loadOWNSProblem
%   cfg.threshold                 see buildEnsembleTransitionThreshold
%   cfg.reduction.perModeTargetRank  reduced rank applied to every mode
%   cfg.reduction.method           per-mode reduction method
%   cfg.angular, cfg.maxDetection, cfg.numerics   as elsewhere
%
% Optional cfg fields:
%   cfg.owns.reducedCacheDirectory  directory for the cached minimal
%                                   {x, A, meta} form of each raw mode
%                                   file, created if missing. Every mode
%                                   is read from its raw file only once;
%                                   all later passes (in this run or a
%                                   later one) reuse the cached copy so
%                                   the expensive high-omega raw files
%                                   are never re-read.
%
% MATLAB version: R2020b

maxSobolDimension = 1111;

fileList = discoverOWNSEnsembleFiles(cfg.owns.ensembleDirectory);
M = numel(fileList);

fprintf('Ensemble first-transition workflow: %d mode files found.\n', M);

% -------------------------------------------------------------------------
% Reference grid selection.
% -------------------------------------------------------------------------

referenceFile = fileList{1};

if isfield(cfg.owns, 'ensembleReferenceFile') && ...
        ~isempty(cfg.owns.ensembleReferenceFile)

    refSelector = cfg.owns.ensembleReferenceFile;

    if isnumeric(refSelector)
        referenceFile = fileList{refSelector};
    else
        referenceFile = char(refSelector);
    end
end

fprintf('Reference grid taken from: %s\n', referenceFile);

referenceMode = loadOWNSModeCached(referenceFile, cfg);
xRef = referenceMode.x;
clear referenceMode;

% -------------------------------------------------------------------------
% Pass 1: cheap mean/variance energy sweep to build the total threshold.
% -------------------------------------------------------------------------

fprintf('Pass 1 of 2: accumulating total mean energy across modes...\n');

overlapMax = xRef(end);
omegaList = zeros(M, 1);
betaList = zeros(M, 1);
fullRankList = zeros(M, 1);

xRefClipped = xRef;

% Only lightweight scalars are retained per mode here; the large A family
% is discarded immediately so memory stays bounded regardless of M. This
% reads (or builds, on first access) the small cached mode form rather
% than the raw OWNS file.
for i = 1:M
    modeData = loadOWNSModeCached(fileList{i}, cfg);

    overlapMax = min(overlapMax, modeData.x(end));

    omegaList(i) = modeData.meta.omega;
    betaList(i) = modeData.meta.beta;
    fullRankList(i) = size(modeData.A{1}, 1);

    clear modeData;
end

xRefClipped = xRefClipped(xRefClipped <= overlapMax + 1e-9);

if numel(xRefClipped) < 2
    error('runOWNSEnsembleWorkflow:InsufficientOverlap', ...
        'The mode files do not share enough streamwise overlap.');
end

meanEnergyTotal = zeros(numel(xRefClipped), 1);
varianceEnergyTotal = zeros(numel(xRefClipped), 1);
meanEnergyByMode = zeros(numel(xRefClipped), M);
interpolatedAByMode = cell(M, 1);

for i = 1:M
    modeData = loadOWNSModeCached(fileList{i}, cfg);

    Ai = interpMatrixFamilyToGrid( ...
        modeData.x, modeData.A, xRefClipped);

    meanEnergyI = cellfun(@trace, Ai);
    varianceEnergyI = cellfun(@(An) 2 * trace(An * An), Ai);

    meanEnergyTotal = meanEnergyTotal + meanEnergyI;
    varianceEnergyTotal = varianceEnergyTotal + varianceEnergyI;
    meanEnergyByMode(:, i) = meanEnergyI;

    interpolatedAByMode{i} = Ai;
    clear modeData Ai;

    if mod(i, max(floor(M / 10), 1)) == 0 || i == M
        fprintf('  pass 1 completed mode %d of %d\n', i, M);
    end
end

[eThreshTotal, thresholdInfo] = buildEnsembleTransitionThreshold( ...
    xRefClipped, meanEnergyTotal, varianceEnergyTotal, cfg.threshold);

fprintf('Total ensemble threshold = %.6e (method: %s)\n', ...
    eThreshTotal(1), thresholdInfo.method);

energyDiagnostics = computeEnsembleEnergyDiagnostics( ...
    xRefClipped, meanEnergyByMode, omegaList, betaList, ...
    fileList, interpolatedAByMode, cfg);

% -------------------------------------------------------------------------
% Pass 2: build and reduce G_i(x) per mode, then assemble block-diagonal.
% -------------------------------------------------------------------------

fprintf('Pass 2 of 2: normalizing and reducing each mode...\n');

perModeCfg = cfg.reduction;
perModeCfg.targetRank = cfg.reduction.perModeTargetRank;
if isfield(perModeCfg, 'verbose')
    perModeCfg.verbose = false;
end

GredByMode = cell(M, 1);
reducedRankList = zeros(M, 1);
reductionInfoByMode = cell(M, 1);

for i = 1:M
    Ai = interpolatedAByMode{i};

    Gi = cell(numel(Ai), 1);
    for n = 1:numel(Ai)
        Gi{n} = Ai{n} / eThreshTotal(n);
    end

    [Gred, ~, redInfo] = reduceEnsembleMode(Gi, xRefClipped, perModeCfg);

    GredByMode{i} = Gred;
    reducedRankList(i) = size(Gred{1}, 1);
    reductionInfoByMode{i} = redInfo;

    interpolatedAByMode{i} = [];
    clear Ai Gi Gred;

    if mod(i, max(floor(M / 10), 1)) == 0 || i == M
        fprintf('  pass 2 completed mode %d of %d\n', i, M);
    end
end

rFull = sum(fullRankList);
rTotal = sum(reducedRankList);

fprintf('Ensemble stochastic dimension: r = %d (from full rank %d)\n', ...
    rTotal, rFull);

if strcmpi(cfg.angular.method, 'rqmc') && rTotal > maxSobolDimension
    error('runOWNSEnsembleWorkflow:SobolDimensionExceeded', ...
        ['Total reduced stochastic dimension r = %d exceeds the Sobol'' ', ...
        'construction limit of %d dimensions. Reduce ', ...
        'cfg.reduction.perModeTargetRank or use cfg.angular.method = ', ...
        '''random''.'], rTotal, maxSobolDimension);
end

[G, blockSizes] = assembleBlockDiagonalFamily(GredByMode);
clear GredByMode;

% -------------------------------------------------------------------------
% Streamwise derivative and probability calculation, reusing the existing
% single-mode machinery unchanged.
% -------------------------------------------------------------------------

meta = struct();
meta.omega = omegaList;
meta.beta = betaList;

Gprime = [];
if ismember(lower(cfg.maxDetection.method), {'hermite', 'adaptive'})
    fprintf('Building Gprime(x)...\n');
    Gprime = buildGprimeFamily(xRefClipped, G, meta, cfg, []);
end

fprintf('Computing angular first-transition integral...\n');

if strcmpi(cfg.angular.method, 'rqmc') && cfg.angular.numReplicates > 1
    prob = computeRQMCCDF( ...
        xRefClipped, G, Gprime, rTotal, ...
        cfg.angular, cfg.maxDetection, cfg.numerics);
else
    [U, wAng, ~] = generateAngularRule(rTotal, cfg.angular);
    validateAngularRule(U, wAng);

    switch lower(cfg.maxDetection.method)
        case 'grid'
            maxData = computeRunningMaxGrid(G, U, cfg);
        case 'hermite'
            maxData = computeRunningMaxHermite(xRefClipped, G, Gprime, U, cfg);
        otherwise
            error('runOWNSEnsembleWorkflow:UnknownMaxDetection', ...
                'Unknown max detection method: %s', cfg.maxDetection.method);
    end

    validateRunningMax(maxData.m);
    prob = computeTransitionCDF(maxData.m, wAng, rTotal, cfg);
end

fprintf('Terminal transition probability F(xmax) = %.6e\n', prob.F(end));

results = struct();
results.x = xRefClipped;
results.prob = prob;
results.meta = meta;
results.thresholdInfo = thresholdInfo;
results.reductionInfoByMode = reductionInfoByMode;
results.blockSizes = blockSizes;
results.rTotal = rTotal;
results.rFull = rFull;
results.numModes = M;
results.fileList = fileList;
results.cfg = cfg;
results.meanEnergyByMode = meanEnergyByMode;
results.omega = omegaList;
results.beta = betaList;
results.energyDiagnostics = energyDiagnostics;

end
