% runOWNSThresholdStudy.m
%
% Explore constant transition thresholds and inlet-amplitude scaling using
% one reusable RQMC energy envelope.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS threshold and amplitude study\n');
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
cfg.angular.numDirections = 4096;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 271828;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

cfg.reduction = struct();
cfg.reduction.method = 'integrated';
cfg.reduction.targetRank = 30;
cfg.reduction.integrationWeight = 'uniform';
cfg.reduction.verbose = false;

cfg.thresholdSweep = struct();
cfg.thresholdSweep.values = [];
cfg.thresholdSweep.factors = logspace(-0.5, 0.5, 9);
cfg.thresholdSweep.amplitudes = [0.5, 0.75, 1.0, 1.25, 1.5];
cfg.thresholdSweep.referenceMethod = 'maximumMeanEnergy';

%% ========================================================================
% Load OWNS data and construct rank-30 integrated basis
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
Afull = problem.A;

% The integrated basis is invariant to a scalar constant threshold.
GforBasis = Afull;

[Vred, reductionInfo] = ...
    buildTransitionBasis( ...
        GforBasis, x, cfg);

Ared = applyReductionToMatrixFamily( ...
    Afull, Vred);

r = size(Vred, 2);

[Aprime, derivativeInfo] = ...
    buildGprimeFamily( ...
        x, Ared, struct(), cfg, []);

fprintf('Reduced rank: %d\n', r);
fprintf('Derivative method: %s\n', ...
    derivativeInfo.method);

%% ========================================================================
% Build reusable energy envelope
% ========================================================================

envelope = buildRQMCEnergyEnvelope( ...
    x, Ared, Aprime, r, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

%% ========================================================================
% Threshold and amplitude family
% ========================================================================

[thresholdValues, thresholdInfo] = ...
    buildThresholdSweepValues( ...
        Ared, cfg.thresholdSweep);

amplitudeValues = ...
    cfg.thresholdSweep.amplitudes;

family = evaluateThresholdAmplitudeFamily( ...
    envelope, ...
    thresholdValues, ...
    amplitudeValues);

fprintf('\nTerminal transition probabilities:\n');

fprintf('%15s', 'threshold');

for j = 1:numel(amplitudeValues)
    fprintf(' %12s', ...
        sprintf('eps=%.2f', amplitudeValues(j)));
end

fprintf('\n');

for i = 1:numel(thresholdValues)
    fprintf('%15.6e', thresholdValues(i));

    for j = 1:numel(amplitudeValues)
        fprintf(' %12.6f', ...
            family.terminalProbability(i, j));
    end

    fprintf('\n');
end

%% ========================================================================
% Solve threshold for target terminal probability
% ========================================================================

targetOptions = struct();
targetOptions.targetIndex = numel(x);
targetOptions.initialThreshold = ...
    thresholdInfo.referenceValue;
targetOptions.probabilityTolerance = 1e-6;
targetOptions.relativeTolerance = 1e-6;
targetOptions.maxIterations = 100;
targetOptions.verbose = true;

targetSolution = ...
    solveThresholdForTargetProbability( ...
        envelope, ...
        0.5, ...
        1.0, ...
        targetOptions);

%% ========================================================================
% Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Threshold and amplitude sensitivity');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
hold on;

colors = lines(numel(amplitudeValues));

for j = 1:numel(amplitudeValues)
    semilogx( ...
        thresholdValues, ...
        family.terminalProbability(:, j), ...
        'o-', ...
        'LineWidth', 1.5, ...
        'Color', colors(j, :), ...
        'DisplayName', ...
        sprintf('\\epsilon = %.2f', ...
            amplitudeValues(j)));
end

xlabel('Constant energy threshold');
ylabel('F_{X_{tr}}(x_{max})');
title('Terminal transition probability');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

for j = 1:numel(amplitudeValues)
    semilogx( ...
        thresholdValues, ...
        family.pCensored(:, j), ...
        'o-', ...
        'LineWidth', 1.5, ...
        'Color', colors(j, :), ...
        'DisplayName', ...
        sprintf('\\epsilon = %.2f', ...
            amplitudeValues(j)));
end

xlabel('Constant energy threshold');
ylabel('S(x_{max})');
title('Censoring probability');
grid on;
box on;

nexttile;
hold on;

% Plot the full CDF at amplitude one for every threshold.
[~, amplitudeOneIndex] = min(abs(amplitudeValues - 1));

thresholdColors = parula(numel(thresholdValues));

for i = 1:numel(thresholdValues)
    plot( ...
        x, ...
        family.F(:, i, amplitudeOneIndex), ...
        'LineWidth', 1.1, ...
        'Color', thresholdColors(i, :), ...
        'DisplayName', ...
        sprintf('e_{th}=%.3e', ...
            thresholdValues(i)));
end

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('CDF sensitivity to threshold');
grid on;
box on;

nexttile;
hold on;

for i = 1:numel(thresholdValues)
    plot( ...
        x, ...
        family.memoryCorrection(:, i, amplitudeOneIndex), ...
        'LineWidth', 1.1, ...
        'Color', thresholdColors(i, :));
end

xlabel('Wall arc length');
ylabel('F_{X_{tr}}-p_{local}');
title('Memory correction');
grid on;
box on;

sgtitle('OWNS threshold and inlet-amplitude sensitivity');

fprintf('\n============================================================\n');
fprintf('OWNS threshold study completed.\n');
fprintf('============================================================\n');