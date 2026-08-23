% runOWNSReductionStudy.m
%
% Compare transition-oriented reduced ranks for the OWNS dataset.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS transition-oriented rank study\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  Configuration
% ========================================================================

cfg = struct();

cfg.gaussianConvention = 'real';

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
cfg.angular.numDirections = 1024;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 314159;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

cfg.reduction = struct();
cfg.reduction.method = 'integrated';
cfg.reduction.targetRank = [];
cfg.reduction.numIntegratedModes = 10;
cfg.reduction.numModesPerStation = 3;
cfg.reduction.stationSelection = 'recordsAndUniform';
cfg.reduction.numUniformStations = 12;
cfg.reduction.recordRelativeTolerance = 1e-6;
cfg.reduction.integrationWeight = 'uniform';
cfg.reduction.useGreedyEnrichment = true;
cfg.reduction.greedyTolerance = 1e-8;
cfg.reduction.maxGreedyModes = 50;
cfg.reduction.snapshotTolerance = 1e-12;
cfg.reduction.verbose = false;

rankValues = [5, 10, 15, 20, 30, 40, 60, 80, 100, 141];

%% ========================================================================
%  Load OWNS data once
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
eThresh = problem.eThresh;
meta = problem.meta;

Gfull = buildGFromAFamily( ...
    A, eThresh, cfg);

rFull = size(Gfull{1}, 1);

%% ========================================================================
%  Construct one nested integrated basis
% ========================================================================

[Kint, integratedInfo] = ...
    buildIntegratedTransitionMatrix( ...
        Gfull, x, cfg.reduction);

[Vall, Dall] = eig(Kint);
lambdaAll = real(diag(Dall));

[lambdaAll, order] = sort( ...
    lambdaAll, 'descend');

Vall = real(Vall(:, order));

[Vall, ~] = qr(Vall, 0);

%% ========================================================================
%  Rank study
% ========================================================================

numRanks = numel(rankValues);

terminalF = zeros(numRanks, 1);
terminalSE = zeros(numRanks, 1);
maximumCdfChange = nan(numRanks, 1);
relativeResidual = zeros(numRanks, 1);
relativeOmittedResidual = zeros(numRanks, 1);
relativeCrossResidual = zeros(numRanks, 1);

allCdf = cell(numRanks, 1);

previousF = [];

fprintf('%8s %16s %14s %14s %14s\n', ...
    'rank', 'F(xmax)', 'SE', 'CDF change', 'G residual');

for j = 1:numRanks
    s = rankValues(j);

    if s > rFull
        continue;
    end

    Vred = Vall(:, 1:s);

    Ared = applyReductionToMatrixFamily( ...
        A, Vred);

    Gred = buildGFromAFamily( ...
        Ared, eThresh, cfg);

    [GprimeRed, ~] = buildGprimeFamily( ...
        x, Gred, struct(), cfg, []);

    prob = computeRQMCCDF( ...
        x, Gred, GprimeRed, s, ...
        cfg.angular, ...
        cfg.maxDetection, ...
        cfg.numerics);

    reductionDiagnostics = ...
        computeReductionDiagnostics( ...
            Gfull, Vred);

    terminalF(j) = prob.F(end);
    terminalSE(j) = prob.terminalStandardError;

    relativeResidual(j) = ...
        reductionDiagnostics.relativeMaximumFullResidual;

    relativeOmittedResidual(j) = ...
        reductionDiagnostics.relativeMaximumOmittedResidual;

    relativeCrossResidual(j) = ...
        reductionDiagnostics.relativeMaximumCrossResidual;

    allCdf{j} = prob.F;

    if ~isempty(previousF)
        maximumCdfChange(j) = ...
            max(abs(prob.F - previousF));
    end

    fprintf('%8d %16.8e %14.6e %14.6e %14.6e\n', ...
        s, ...
        terminalF(j), ...
        terminalSE(j), ...
        maximumCdfChange(j), ...
        relativeResidual(j));

    previousF = prob.F;
end

%% ========================================================================
%  Integrated spectrum
% ========================================================================

cumulativeIntegratedFraction = ...
    cumsum(max(lambdaAll, 0)) ...
    / sum(max(lambdaAll, 0));

fprintf('\nIntegrated-matrix cumulative fractions:\n');

for s = rankValues
    if s <= rFull
        fprintf('  rank %3d: %.8f\n', ...
            s, cumulativeIntegratedFraction(s));
    end
end

%% ========================================================================
%  Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'OWNS stochastic-rank convergence');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

errorbar( ...
    rankValues, terminalF, ...
    2 * terminalSE, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Reduced stochastic rank');
ylabel('F_{X_{tr}}(x_{max})');
title('Terminal transition probability');
grid on;
box on;

nexttile;

semilogy( ...
    rankValues, relativeResidual, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

hold on;

semilogy( ...
    rankValues, relativeOmittedResidual, ...
    's--', ...
    'LineWidth', 1.25);

semilogy( ...
    rankValues, relativeCrossResidual, ...
    'd-.', ...
    'LineWidth', 1.25);

xlabel('Reduced stochastic rank');
ylabel('Relative matrix residual');
title('Projection residuals');
legend('Full', 'Omitted block', 'Cross block', ...
    'Location', 'best');
grid on;
box on;

nexttile;

semilogy( ...
    rankValues(2:end), ...
    maximumCdfChange(2:end), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Reduced stochastic rank');
ylabel('Maximum CDF change');
title('Rank-to-rank CDF convergence');
grid on;
box on;

nexttile;

semilogy( ...
    1:rFull, ...
    max(lambdaAll, eps), ...
    'LineWidth', 1.5);

xlabel('Integrated mode index');
ylabel('Integrated eigenvalue');
title('Integrated transition spectrum');
grid on;
box on;

sgtitle('OWNS transition-oriented dimension reduction');

figure('Color', 'w', ...
    'Name', 'OWNS CDF rank comparison');

hold on;

colors = lines(numRanks);

for j = 1:numRanks
    if ~isempty(allCdf{j})
        plot(x, allCdf{j}, ...
            'LineWidth', 1.25, ...
            'Color', colors(j, :), ...
            'DisplayName', ...
            sprintf('s = %d', rankValues(j)));
    end
end

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('First-transition CDF by reduced rank');
legend('Location', 'eastoutside');
grid on;
box on;

fprintf('\n============================================================\n');
fprintf('OWNS rank study completed.\n');
fprintf('============================================================\n');
