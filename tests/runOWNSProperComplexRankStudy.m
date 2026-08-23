% runOWNSProperComplexRankStudy.m
%
% Verify proper-complex rank convergence using a complex-specific
% integrated basis.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS proper-complex rank study\n');
fprintf('============================================================\n\n');

%% Configuration

cfg = struct();

cfg.owns = struct();
cfg.owns.dataFile = ...
    '/data2/jampauls/DataforFigures/SaveData/FirstPassingCode/data/ninth_Run_i_1_j_1_solution.mat';
cfg.owns.solutionVariable = 'solution';
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;
cfg.owns.weightTolerance = 1e-12;
cfg.owns.verbose = false;

cfg.coordinates = struct();
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.referenceIndex = 1;
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;

cfg.numerics = struct();
cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.maxDetection = struct();
cfg.maxDetection.method = 'hermite';

cfg.angular = struct();
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 2048;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 424242;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

basisOptions = struct();
basisOptions.integrationWeight = 'uniform';
basisOptions.verbose = true;

rankValues = [10, 20, 30, 40, 60, 80, 100, 141];

%% Load once and construct both A families

loaded = load( ...
    cfg.owns.dataFile, ...
    cfg.owns.solutionVariable);

solution = loaded.(cfg.owns.solutionVariable);
clear loaded;

[x, coordinateMeta] = ...
    extractOWNSCoordinate( ...
        solution, cfg.coordinates);

[Areal, matrixInfo, Acomplex] = ...
    buildAFromOWNS(solution, cfg);

clear solution;

threshold = max(matrixInfo.meanEnergy);
eThresh = repmat(threshold, numel(x), 1);

Gcomplex = buildGFromAFamily( ...
    Acomplex, eThresh, cfg);

validateGFamily(Gcomplex, cfg);

[GprimeComplex, derivativeInfo] = ...
    buildGprimeFamily( ...
        x, Gcomplex, struct(), cfg, []);

validateGprimeFamily( ...
    x, Gcomplex, GprimeComplex);

fprintf('Derivative method: %s\n', ...
    derivativeInfo.method);

%% Full ordered complex integrated basis

[Vall, basisInfo] = ...
    buildComplexIntegratedBasis( ...
        Gcomplex, x, size(Gcomplex{1}, 1), ...
        basisOptions);

%% Coupled rank calculation

study = computeCoupledComplexRankStudy( ...
    x, Gcomplex, GprimeComplex, ...
    Vall, rankValues, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

%% Print summary

fprintf('\nProper-complex rank summary:\n');
fprintf('%8s %16s %14s %16s %14s\n', ...
    'rank', ...
    'F(xmax)', ...
    'SE', ...
    'diff full', ...
    'paired SE');

for j = 1:numel(rankValues)
    fprintf('%8d %16.8e %14.6e %16.6e %14.6e\n', ...
        rankValues(j), ...
        study.FMean(end, j), ...
        study.FStandardError(end, j), ...
        study.maximumReferenceDifference(j), ...
        study.maximumReferenceDifferenceSE(j));
end

fprintf('\nIntegrated complex spectrum:\n');

for j = 1:numel(rankValues)
    s = rankValues(j);

    fprintf('  rank %3d: %.10f\n', ...
        s, basisInfo.cumulativeFraction(s));
end

%% Compare shared-real and complex-specific rank-30 bases

realBasisCfg = struct();
realBasisCfg.reduction = struct();
realBasisCfg.reduction.method = 'integrated';
realBasisCfg.reduction.targetRank = 30;
realBasisCfg.reduction.integrationWeight = 'uniform';
realBasisCfg.reduction.verbose = false;

Greal = buildGFromAFamily( ...
    Areal, eThresh, cfg);

[Vreal30, ~] = buildTransitionBasis( ...
    Greal, x, realBasisCfg);

Vcomplex30 = Vall(:, 1:30);

GcomplexRealBasis = ...
    applyComplexReductionToMatrixFamily( ...
        Gcomplex, Vreal30);

GprimeComplexRealBasis = ...
    applyComplexReductionToMatrixFamily( ...
        GprimeComplex, Vreal30);

GcomplexOwnBasis = ...
    applyComplexReductionToMatrixFamily( ...
        Gcomplex, Vcomplex30);

GprimeComplexOwnBasis = ...
    applyComplexReductionToMatrixFamily( ...
        GprimeComplex, Vcomplex30);

resultSharedRealBasis = ...
    computeProperComplexRQMCCDF( ...
        x, ...
        GcomplexRealBasis, ...
        GprimeComplexRealBasis, ...
        30, ...
        cfg.angular, ...
        cfg.maxDetection, ...
        cfg.numerics);

resultComplexBasis = ...
    computeProperComplexRQMCCDF( ...
        x, ...
        GcomplexOwnBasis, ...
        GprimeComplexOwnBasis, ...
        30, ...
        cfg.angular, ...
        cfg.maxDetection, ...
        cfg.numerics);

fprintf('\nRank-30 basis comparison:\n');
fprintf('  shared real basis F(xmax) = %.12e +/- %.3e\n', ...
    resultSharedRealBasis.F(end), ...
    2 * resultSharedRealBasis.terminalStandardError);

fprintf('  complex basis F(xmax)     = %.12e +/- %.3e\n', ...
    resultComplexBasis.F(end), ...
    2 * resultComplexBasis.terminalStandardError);

fprintf('  full complex F(xmax)      = %.12e +/- %.3e\n', ...
    study.FMean(end, end), ...
    2 * study.FStandardError(end, end));

%% Plots

figure('Color', 'w', ...
    'Name', 'Proper-complex rank convergence');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

errorbar( ...
    rankValues, ...
    study.FMean(end, :), ...
    2 * study.FStandardError(end, :), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Complex reduced rank');
ylabel('F_{X_{tr}}(x_{max})');
title('Terminal transition probability');
grid on;
box on;

nexttile;

referenceY = ...
    study.maximumReferenceDifference(1:end - 1);

referenceSE = ...
    study.maximumReferenceDifferenceSE(1:end - 1);

referenceRanks = rankValues(1:end - 1);

lowerError = min( ...
    2 * referenceSE, ...
    0.95 * referenceY);

upperError = 2 * referenceSE;

valid = referenceY > 0 & isfinite(referenceY);

errorbar( ...
    referenceRanks(valid), ...
    referenceY(valid), ...
    lowerError(valid), ...
    upperError(valid), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

set(gca, 'YScale', 'log');

xlabel('Complex reduced rank');
ylabel('Maximum paired CDF difference');
title('Difference from full complex rank');
grid on;
box on;

nexttile;

semilogy( ...
    1:numel(basisInfo.eigenvalues), ...
    max(basisInfo.eigenvalues, eps), ...
    'LineWidth', 1.5);

xlabel('Mode index');
ylabel('Integrated eigenvalue');
title('Complex integrated spectrum');
grid on;
box on;

nexttile;
hold on;

plot(x, resultSharedRealBasis.F, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Rank 30, real basis');

plot(x, resultComplexBasis.F, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Rank 30, complex basis');

plot(x, study.FMean(:, end), ':', ...
    'LineWidth', 2, ...
    'DisplayName', 'Full complex rank');

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('Basis sensitivity');
legend('Location', 'best');
grid on;
box on;

sgtitle('Proper-complex OWNS verification');

fprintf('\n============================================================\n');
fprintf('Proper-complex rank study completed.\n');
fprintf('============================================================\n');