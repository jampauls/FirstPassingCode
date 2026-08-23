% runOWNSGaussianConventionStudy.m
%
% Compare real Gaussian coefficients with proper-complex coefficients for
% the same OWNS propagated factor and physical energy threshold.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS Gaussian-coordinate convention study\n');
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
cfg.owns.verbose = true;

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
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;

cfg.angular = struct();
cfg.angular.method = 'rqmc';
cfg.angular.numDirections = 4096;
cfg.angular.numReplicates = 8;
cfg.angular.randomSeed = 161803;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

cfg.reduction = struct();
cfg.reduction.method = 'integrated';
cfg.reduction.targetRank = 30;
cfg.reduction.integrationWeight = 'uniform';
cfg.reduction.verbose = false;

%% ========================================================================
% Load the OWNS solution once
% ========================================================================

loaded = load( ...
    cfg.owns.dataFile, ...
    cfg.owns.solutionVariable);

solution = ...
    loaded.(cfg.owns.solutionVariable);

clear loaded;

inspectOWNSData(solution, cfg);

[x, coordinateMeta] = ...
    extractOWNSCoordinate( ...
        solution, cfg.coordinates);

[ArealFull, matrixInfo, AcomplexFull] = ...
    buildAFromOWNS(solution, cfg);

clear solution;

%% ========================================================================
% Use the same physical constant threshold
% ========================================================================

% The mean energy is the same in both coefficient conventions.
threshold = max(matrixInfo.meanEnergy);

eThresh = repmat(threshold, numel(x), 1);

fprintf('Comparison threshold = %.12e\n\n', threshold);

%% ========================================================================
% Construct one shared real integrated basis
% ========================================================================

% Use the current real-coefficient model to construct the fixed rank-30
% basis. Apply that same basis to both stochastic models so the comparison
% isolates the coefficient convention.
GrealFull = buildGFromAFamily( ...
    ArealFull, eThresh, cfg);

[Vred, reductionInfo] = ...
    buildTransitionBasis( ...
        GrealFull, x, cfg);

Areal = applyReductionToMatrixFamily( ...
    ArealFull, Vred);

Acomplex = applyReductionToMatrixFamily( ...
    AcomplexFull, Vred);

r = size(Vred, 2);

Greal = buildGFromAFamily( ...
    Areal, eThresh, cfg);

Gcomplex = buildGFromAFamily( ...
    Acomplex, eThresh, cfg);

[GprimeReal, ~] = buildGprimeFamily( ...
    x, Greal, struct(), cfg, []);

[GprimeComplex, ~] = buildGprimeFamily( ...
    x, Gcomplex, struct(), cfg, []);

%% ========================================================================
% Real-coefficient probability
% ========================================================================

realResult = computeRQMCCDF( ...
    x, Greal, GprimeReal, r, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

%% ========================================================================
% Proper-complex coefficient probability
% ========================================================================

complexResult = computeProperComplexRQMCCDF( ...
    x, Gcomplex, GprimeComplex, r, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

%% ========================================================================
% Moment comparison
% ========================================================================

meanReal = zeros(numel(x), 1);
meanComplex = zeros(numel(x), 1);

varianceReal = zeros(numel(x), 1);
varianceComplex = zeros(numel(x), 1);

for n = 1:numel(x)
    meanReal(n) = trace(Areal{n});
    meanComplex(n) = real(trace(Acomplex{n}));

    varianceReal(n) = ...
        2 * trace(Areal{n} * Areal{n});

    varianceComplex(n) = ...
        real(trace(Acomplex{n} * Acomplex{n}));
end

meanDifference = max(abs( ...
    meanReal - meanComplex));

fprintf('Gaussian convention comparison:\n');
fprintf('  reduced rank               = %d\n', r);
fprintf('  maximum mean difference    = %.12e\n', ...
    meanDifference);
fprintf('  real F(xmax)               = %.12e\n', ...
    realResult.F(end));
fprintf('  real terminal SE           = %.12e\n', ...
    realResult.terminalStandardError);
fprintf('  proper-complex F(xmax)     = %.12e\n', ...
    complexResult.F(end));
fprintf('  proper-complex terminal SE = %.12e\n', ...
    complexResult.terminalStandardError);
fprintf('  terminal probability diff  = %.12e\n\n', ...
    complexResult.F(end) - realResult.F(end));

%% ========================================================================
% Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Gaussian-coordinate convention comparison');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
hold on;

plot(x, realResult.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'Real coefficients');

plot(x, complexResult.F, '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Proper-complex coefficients');

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('First-transition CDF');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

plot(x, realResult.pLocal, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Real coefficients');

plot(x, complexResult.pLocal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Proper-complex coefficients');

xlabel('Wall arc length');
ylabel('Local exceedance probability');
title('Local energy tails');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

plot(x, varianceReal, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Real coefficient variance');

plot(x, varianceComplex, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Proper-complex variance');

xlabel('Wall arc length');
ylabel('Energy variance');
title('One-point energy variance');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, ...
    complexResult.F - realResult.F, ...
    'LineWidth', 1.5);

xlabel('Wall arc length');
ylabel('\Delta F');
title('Proper-complex minus real');
grid on;
box on;

sgtitle('Sensitivity to Gaussian coefficient convention');

fprintf('============================================================\n');
fprintf('Gaussian convention study completed.\n');
fprintf('============================================================\n');
