% runOWNSPositiveIncrementStudy.m
%
% Diagnose whether significant positive energy increments occupy a small
% fixed stochastic subspace, leaving a monotonically decaying complement.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS positive-increment growth-subspace study\n');
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

cfg.numerics = struct();
cfg.numerics.symmetrizeG = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.growthSubspace = struct();

cfg.growthSubspace.incrementScaling = 'perUnitX';

cfg.growthSubspace.positiveEigenvalueTolerances = ...
    [1e-6, 1e-8, 1e-10, 1e-12, 1e-14];

cfg.growthSubspace.positivePartTolerance = 1e-12;

cfg.growthSubspace.aggregateWeightMethod = 'magnitude';

cfg.growthSubspace.candidateRanks = ...
    [1, 2, 3, 5, 10, 15, 20, 30, 40, 60, 80, 100, 120, 141];

cfg.growthSubspace.complementTolerance = 1e-10;
cfg.growthSubspace.verbose = true;

% Optional directional record-growth diagnostic.
cfg.recordGrowth = struct();
cfg.recordGrowth.enable = true;
cfg.recordGrowth.numDirections = 8192;
cfg.recordGrowth.randomSeed = 8675309;

% -------------------------------------------------------------------------
% Threshold and inlet-amplitude exploration
% -------------------------------------------------------------------------

cfg.thresholdSweep = struct();

cfg.thresholdSweep.enable = false;

% Constant threshold values. If empty, values can be generated relative
% to the maximum mean energy.
cfg.thresholdSweep.values = [];

% Multiplicative factors applied to a reference threshold.
cfg.thresholdSweep.factors = ...
    logspace(-0.5, 0.5,  9);

% Inlet-factor amplitude multipliers:
%
%   B0(epsilon) = epsilon * B0(reference)
%
% Energy therefore scales by epsilon^2.
cfg.thresholdSweep.amplitudes = [0.5, 0.75, 1.0, 1.25, 1.5];

% If thresholdSweep.values is empty, use:
%
%   threshold = factor * max_x E[e(x)].
cfg.thresholdSweep.referenceMethod = 'maximumMeanEnergy';

% Optional threshold solution for a target terminal probability.
cfg.thresholdSweep.solveTarget.enable = false;
cfg.thresholdSweep.solveTarget.probability = 0.5;
cfg.thresholdSweep.solveTarget.amplitude = 1.0;
cfg.thresholdSweep.solveTarget.relativeTolerance = 1e-6;
cfg.thresholdSweep.solveTarget.probabilityTolerance = 1e-6;
cfg.thresholdSweep.solveTarget.maxIterations = 100;


%% ========================================================================
% Load and build A(x)
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
meta = problem.meta;

clear problem;

%% ========================================================================
% Positive-increment analysis
% ========================================================================

growthReport = analyzePositiveIncrementSubspace( ...
    A, x, cfg.growthSubspace);

%% ========================================================================
% Optional sampled historical-record growth
% ========================================================================

recordReport = [];

if cfg.recordGrowth.enable
    angularCfg = struct();

    angularCfg.method = 'rqmc';
    angularCfg.numDirections = ...
        cfg.recordGrowth.numDirections;

    angularCfg.randomSeed = ...
        cfg.recordGrowth.randomSeed;

    angularCfg.replicateIndex = 1;
    angularCfg.sobolSkip = 1024;
    angularCfg.sobolLeap = 0;
    angularCfg.scramble = true;
    angularCfg.probabilityClip = 1e-12;

    r = size(A{1}, 1);

    [U, weights] = generateAngularRule( ...
        r, angularCfg);

    recordReport = analyzeDirectionalRecordGrowth( ...
        A, U, weights, x);
end

%% ========================================================================
% Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Positive increment structure');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% Number of positive directions
% -------------------------------------------------------------------------

nexttile;
hold on;

colors = lines(numel( ...
    growthReport.positiveEigenvalueTolerances));

for iTol = 1:numel( ...
        growthReport.positiveEigenvalueTolerances)

    plot( ...
        growthReport.intervalCoordinate, ...
        growthReport.positiveRank(:, iTol), ...
        'LineWidth', 1.25, ...
        'Color', colors(iTol, :), ...
        'DisplayName', ...
        sprintf('\\tau=%.0e', ...
            growthReport.positiveEigenvalueTolerances(iTol)));
end

xlabel('Wall arc length');
ylabel('Positive increment rank');
title('Locally growing stochastic directions');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Positive and negative variation
% -------------------------------------------------------------------------

nexttile;
hold on;

semilogy( ...
    growthReport.intervalCoordinate, ...
    max(growthReport.positiveTrace, realmin), ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Positive trace');

semilogy( ...
    growthReport.intervalCoordinate, ...
    max(growthReport.negativeTrace, realmin), ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Negative trace');

xlabel('Wall arc length');
ylabel('Increment trace magnitude');
title('Positive and negative energy variation');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Cumulative positive-range dimension
% -------------------------------------------------------------------------

nexttile;
hold on;

for iTol = 1:numel( ...
        growthReport.positiveEigenvalueTolerances)

    plot( ...
        growthReport.intervalCoordinate, ...
        growthReport.cumulativePositiveRank(:, iTol), ...
        'LineWidth', 1.5, ...
        'Color', colors(iTol, :), ...
        'DisplayName', ...
        sprintf('\\tau=%.0e', ...
            growthReport.positiveEigenvalueTolerances(iTol)));
end

xlabel('Wall arc length');
ylabel('Cumulative growth-subspace rank');
title('Historical span of positive increments');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Aggregate growth spectrum
% -------------------------------------------------------------------------

nexttile;

semilogy( ...
    1:numel(growthReport.growthEigenvalues), ...
    max(growthReport.growthEigenvalues, realmin), ...
    'o-', ...
    'LineWidth', 1.25, ...
    'MarkerSize', 3);

xlabel('Growth mode index');
ylabel('Aggregate positive-growth eigenvalue');
title('Fixed growth-basis spectrum');
grid on;
box on;

sgtitle('Positive energy-increment diagnostics');

figure('Color', 'w', ...
    'Name', 'Growth-subspace quality');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% Aggregate capture
% -------------------------------------------------------------------------

nexttile;
hold on;

plot( ...
    growthReport.candidateRanks, ...
    growthReport.aggregateCapturedFraction, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Aggregate spectrum');

plot( ...
    growthReport.candidateRanks, ...
    growthReport.weightedPositiveTraceCapture, ...
    's--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Weighted positive trace');

xlabel('Fixed growth-subspace dimension');
ylabel('Captured fraction');
title('Positive-growth capture');
legend('Location', 'best');
grid on;
box on;
ylim([0, 1.01]);

% -------------------------------------------------------------------------
% Positive-part reconstruction residual
% -------------------------------------------------------------------------

nexttile;
hold on;

semilogy( ...
    growthReport.candidateRanks, ...
    growthReport.maximumPositiveOperatorResidual, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Maximum');

semilogy( ...
    growthReport.candidateRanks, ...
    growthReport.medianPositiveOperatorResidual, ...
    's--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Median');

xlabel('Fixed growth-subspace dimension');
ylabel('Relative positive-part residual');
title('Positive-increment reconstruction');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Complement growth
% -------------------------------------------------------------------------

nexttile;

semilogy( ...
    growthReport.candidateRanks, ...
    max( ...
        growthReport.maximumRelativeComplementGrowth, ...
        realmin), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Fixed growth-subspace dimension');
ylabel('Maximum relative complement growth');
title('Does the complement still grow?');
grid on;
box on;

% -------------------------------------------------------------------------
% Active-complement coupling
% -------------------------------------------------------------------------

nexttile;

semilogy( ...
    growthReport.candidateRanks, ...
    max( ...
        growthReport.maximumRelativeCrossCoupling, ...
        realmin), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Fixed growth-subspace dimension');
ylabel('Maximum relative cross coupling');
title('Growth/complement coupling');
grid on;
box on;

sgtitle('Quality of fixed growth-coordinate models');

% -------------------------------------------------------------------------
% Sampled directional record-growth plots
% -------------------------------------------------------------------------

if ~isempty(recordReport)
    figure('Color', 'w', ...
        'Name', 'Directional record growth');

    tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    nexttile;

    plot( ...
        x, ...
        recordReport.weightedRecordIncrementByStation, ...
        'LineWidth', 1.5);

    xlabel('Wall arc length');
    ylabel('Weighted record increment');
    title('Sampled historical-envelope growth');
    grid on;
    box on;

    nexttile;

    semilogy( ...
        1:numel(recordReport.recordEigenvalues), ...
        max(recordReport.recordEigenvalues, realmin), ...
        'o-', ...
        'LineWidth', 1.25, ...
        'MarkerSize', 3);

    xlabel('Record-growth mode index');
    ylabel('Sampled record-growth eigenvalue');
    title('Input directions associated with new records');
    grid on;
    box on;
end

fprintf('\n============================================================\n');
fprintf('OWNS positive-increment study completed.\n');
fprintf('============================================================\n');
