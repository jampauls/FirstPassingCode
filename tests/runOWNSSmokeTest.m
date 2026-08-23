% runOWNSSmokeTest.m
%
% Initial end-to-end OWNS data test.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS first-transition smoke test\n');
fprintf('============================================================\n\n');

cfg = struct();

cfg.gaussianConvention = 'real';

cfg.owns = struct();
cfg.owns.dataFile = ...
    '/data2/jampauls/DataforFigures/SaveData/FirstPassingCode/data/ninth_Run_i_1_j_1_solution.mat';
cfg.owns.solutionVariable = 'solution';
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;
cfg.owns.validateWeights = true;
cfg.owns.weightTolerance = 1e-12;
cfg.owns.verbose = true;

cfg.coordinates = struct();
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.referenceIndex = 1;
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;

cfg.threshold = struct();
cfg.threshold.method = 'auto';
cfg.threshold.autoMaxMeanFactor = 1;
cfg.threshold.warnIfProvisional = true;

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
cfg.angular.numDirections = 512;
cfg.angular.numReplicates = 4;
cfg.angular.randomSeed = 1234;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

%% Load and compress

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
eThresh = problem.eThresh;
meta = problem.meta;

fprintf('Compressed matrix family:\n');
fprintf('  number of matrices = %d\n', numel(A));
fprintf('  matrix size        = %d-by-%d\n', ...
    size(A{1}, 1), size(A{1}, 2));
fprintf('  threshold          = %.6e to %.6e\n\n', ...
    min(eThresh), max(eThresh));

%% Construct G and Gprime

G = buildGFromAFamily(A, eThresh, cfg);
validateGFamily(G, cfg);

[Gprime, derivativeInfo] = ...
    buildGprimeFamily(x, G, meta, cfg, []);

validateGprimeFamily(x, G, Gprime);

fprintf('Gprime method: %s\n\n', ...
    derivativeInfo.method);

%% RQMC first-transition calculation

r = size(G{1}, 1);

prob = computeRQMCCDF( ...
    x, G, Gprime, r, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

validateProbabilityCurves(prob, cfg);

fprintf('OWNS smoke-test result:\n');
fprintf('  F(xmax)      = %.12e\n', prob.F(end));
fprintf('  S(xmax)      = %.12e\n', prob.S(end));
fprintf('  terminal SE  = %.12e\n', ...
    prob.terminalStandardError);
fprintf('  max local-CDF violation = %.12e\n\n', ...
    max(prob.pLocal - prob.F));

%% Basic plots

figure('Color', 'w', ...
    'Name', 'OWNS first-transition smoke test');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
plot(x, meta.matrixInfo.meanEnergy, ...
    'LineWidth', 1.5);
hold on;
plot(x, eThresh, '--', ...
    'LineWidth', 1.5);
xlabel('Wall arc length');
ylabel('Energy');
title('Mean energy and provisional threshold');
legend('Mean energy', 'Threshold', ...
    'Location', 'best');
grid on;
box on;

nexttile;
plot(x, prob.F, ...
    'LineWidth', 2);
hold on;
plot(x, prob.pLocal, '--', ...
    'LineWidth', 1.5);
xlabel('Wall arc length');
ylabel('Probability');
title('Local and first-transition probabilities');
legend('First-transition CDF', 'Local exceedance', ...
    'Location', 'best');
grid on;
box on;

nexttile;
plot(x, prob.S, ...
    'LineWidth', 2);
xlabel('Wall arc length');
ylabel('S(x)');
title('Survival probability');
grid on;
box on;

nexttile;
plot(x, prob.memoryCorrection, ...
    'LineWidth', 2);
xlabel('Wall arc length');
ylabel('F - p_{local}');
title('Memory correction');
grid on;
box on;

sgtitle('OWNS first-transition smoke test');

fprintf('============================================================\n');
fprintf('OWNS smoke test completed successfully.\n');
fprintf('============================================================\n');