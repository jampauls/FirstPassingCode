% runOWNSCoupledReductionStudy.m
%
% Coupled RQMC convergence study for nested integrated bases.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('Coupled OWNS stochastic-rank study\n');
fprintf('============================================================\n\n');

%% ========================================================================
% Configuration
% ========================================================================

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

cfg.threshold = struct();
cfg.threshold.method = 'auto';
cfg.threshold.autoMaxMeanFactor = 1;
cfg.threshold.warnIfProvisional = false;

cfg.numerics = struct();
cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.maxDetection = struct();
cfg.maxDetection.method = 'hermite';
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;

cfg.angular = struct();
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 2048;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 314159;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

cfg.reduction = struct();
cfg.reduction.integrationWeight = 'uniform';

rankValues = [5, 10, 15, 20, 30, 40, 60, 80, 100, 141];

%% ========================================================================
% Load and build full matrix family
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
eThresh = problem.eThresh;

Gfull = buildGFromAFamily( ...
    A, eThresh, cfg);

[GprimeFull, derivativeInfo] = ...
    buildGprimeFamily( ...
        x, Gfull, problem.meta, cfg, []);

fprintf('Gprime method: %s\n', derivativeInfo.method);

rFull = size(Gfull{1}, 1);

%% ========================================================================
% Nested integrated basis
% ========================================================================

[Kint, integratedInfo] = ...
    buildIntegratedTransitionMatrix( ...
        Gfull, x, cfg.reduction);

[Vall, Dall] = eig(Kint);

lambdaAll = real(diag(Dall));

[lambdaAll, order] = ...
    sort(lambdaAll, 'descend');

Vall = real(Vall(:, order));

% Eigenvectors of a real symmetric matrix should already be orthonormal.
% Apply only a sign convention, rather than a second factorization that
% could obscure ordering.
for j = 1:size(Vall, 2)
    [~, indexLargest] = max(abs(Vall(:, j)));

    if Vall(indexLargest, j) < 0
        Vall(:, j) = -Vall(:, j);
    end
end

validateReductionBasis(Vall, rFull);

%% ========================================================================
% Coupled study
% ========================================================================

study = computeCoupledRankStudy( ...
    x, Gfull, GprimeFull, Vall, rankValues, ...
    cfg.angular, cfg.maxDetection, cfg.numerics);

%% ========================================================================
% Automatic rank recommendation
% ========================================================================

rankSelectionOptions = struct();

% Desired maximum absolute error in the complete transition CDF.
rankSelectionOptions.cdfTolerance = 5e-3;

% Use a two-standard-error adjustment for the paired RQMC difference.
rankSelectionOptions.confidenceMultiplier = 2;

% Require the normalized matrix-family residual to be below one percent.
rankSelectionOptions.matrixResidualTolerance = 1e-2;

rankSelectionOptions.excludeReferenceRank = true;

rankSelection = selectRankFromCoupledStudy( ...
    study, rankSelectionOptions);

%% ========================================================================
% Print summary
% ========================================================================

fprintf('\nCoupled rank summary:\n');

fprintf([ ...
    '%8s %15s %12s %15s %12s ', ...
    '%15s %12s %13s\n'], ...
    'rank', ...
    'F(xmax)', ...
    'SE', ...
    'diff previous', ...
    'paired SE', ...
    'diff full', ...
    'paired SE', ...
    'G residual');

for j = 1:numel(rankValues)
    fprintf([ ...
        '%8d %15.8e %12.4e %15.4e %12.4e ', ...
        '%15.4e %12.4e %13.4e\n'], ...
        rankValues(j), ...
        study.FMean(end, j), ...
        study.FStandardError(end, j), ...
        study.maximumSuccessiveDifference(j), ...
        study.maximumSuccessiveDifferenceSE(j), ...
        study.maximumReferenceDifference(j), ...
        study.maximumReferenceDifferenceSE(j), ...
        study.relativeFullResidual(j));
end

%% ========================================================================
% Integrated spectrum
% ========================================================================

lambdaPositive = max(lambdaAll, 0);

cumulativeFraction = ...
    cumsum(lambdaPositive) / sum(lambdaPositive);

fprintf('\nIntegrated spectrum:\n');

for j = 1:numel(rankValues)
    s = rankValues(j);

    fprintf('  rank %3d: %.10f\n', ...
        s, cumulativeFraction(s));
end

%% ========================================================================
% Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Coupled OWNS rank convergence');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

successiveX = rankValues(2:end);
successiveY = ...
    study.maximumSuccessiveDifference(2:end);

successiveSE = ...
    study.maximumSuccessiveDifferenceSE(2:end);

% Prevent the lower error-bar endpoint from reaching zero or becoming
% negative on a logarithmic axis.
successiveLowerError = min( ...
    2 * successiveSE, ...
    0.95 * successiveY);

successiveUpperError = ...
    2 * successiveSE;

valid = ...
    isfinite(successiveY) ...
    & isfinite(successiveLowerError) ...
    & isfinite(successiveUpperError) ...
    & successiveY > 0;

errorbar( ...
    successiveX(valid), ...
    successiveY(valid), ...
    successiveLowerError(valid), ...
    successiveUpperError(valid), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

set(gca, 'YScale', 'log');

xlabel('Reduced stochastic rank');
ylabel('Maximum paired CDF difference');
title('Difference from previous rank');
grid on;
box on;


nexttile;

referenceX = rankValues(1:end - 1);
referenceY = ...
    study.maximumReferenceDifference(1:end - 1);

referenceSE = ...
    study.maximumReferenceDifferenceSE(1:end - 1);

referenceLowerError = min( ...
    2 * referenceSE, ...
    0.95 * referenceY);

referenceUpperError = ...
    2 * referenceSE;

valid = ...
    isfinite(referenceY) ...
    & isfinite(referenceLowerError) ...
    & isfinite(referenceUpperError) ...
    & referenceY > 0;

errorbar( ...
    referenceX(valid), ...
    referenceY(valid), ...
    referenceLowerError(valid), ...
    referenceUpperError(valid), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

set(gca, 'YScale', 'log');

xlabel('Reduced stochastic rank');
ylabel('Maximum difference from full rank');
title('Paired reduction error');
grid on;
box on;

nexttile;

semilogy( ...
    rankValues, ...
    study.relativeFullResidual, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

hold on;

semilogy( ...
    rankValues, ...
    study.relativeOmittedResidual, ...
    's--', ...
    'LineWidth', 1.25);

semilogy( ...
    rankValues, ...
    study.relativeCrossResidual, ...
    'd-.', ...
    'LineWidth', 1.25);

xlabel('Reduced stochastic rank');
ylabel('Relative matrix residual');
title('Projection residuals');
legend('Full', 'Omitted', 'Cross', ...
    'Location', 'best');
grid on;
box on;

sgtitle('Coupled transition-oriented rank convergence');

figure('Color', 'w', ...
    'Name', 'Coupled rank CDFs');

hold on;

colors = lines(numel(rankValues));

for j = 1:numel(rankValues)
    plot(x, study.FMean(:, j), ...
        'LineWidth', 1.25, ...
        'Color', colors(j, :), ...
        'DisplayName', ...
        sprintf('s = %d', rankValues(j)));
end

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('Coupled RQMC transition CDFs');
legend('Location', 'eastoutside');
grid on;
box on;

fprintf('\n============================================================\n');
fprintf('Coupled OWNS rank study completed.\n');
fprintf('============================================================\n');