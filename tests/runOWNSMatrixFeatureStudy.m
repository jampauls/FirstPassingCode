% runOWNSMatrixFeatureStudy.m
%
% Analyze intrinsic matrix-space dimension of:
%
%   1. A_k;
%   2. A_k - A_0;
%   3. positive energy increments D_k^+.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS matrix-trajectory feature study\n');
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

% Structural study does not use a physical threshold.
cfg.threshold = struct();
cfg.threshold.method = 'specifiedScalar';
cfg.threshold.value = 1;
cfg.threshold.warnIfProvisional = false;

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
cfg.growthSubspace.verbose = false;

cfg.matrixFeatures = struct();

cfg.matrixFeatures.candidateDimensions = ...
    [0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 60, ...
     80, 100, 120, 140, 160, 180, 200, 220, 240, 260];

cfg.matrixFeatures.rankTolerances = ...
    [1e-6, 1e-8, 1e-10, 1e-12, 1e-14];

cfg.matrixFeatures.spectralErrorTargets = ...
    [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4, 1e-6];

cfg.matrixFeatures.numBasisMatricesToStore = 60;
cfg.matrixFeatures.verbose = true;

cfg.output = struct();

cfg.output.file = fullfile( ...
    projectRoot, ...
    'results', ...
    'owns_matrix_feature_study.mat');

%% ========================================================================
% Load matrix family
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
meta = problem.meta;

clear problem;

Nx = numel(A);
r = size(A{1}, 1);

%% ========================================================================
% Linear family: A_k
% ========================================================================

linearOptions = cfg.matrixFeatures;
linearOptions.name = 'linear family A_k';

linearReport = analyzeMatrixTrajectoryFeatures( ...
    A, linearOptions);

%% ========================================================================
% Affine family: A_k - A_0
% ========================================================================

A0 = 0.5 * (A{1} + A{1}.');

affineFamily = cell(Nx, 1);

for k = 1:Nx
    affineFamily{k} = ...
        0.5 * (A{k} + A{k}.') ...
        - A0;

    affineFamily{k} = ...
        0.5 * ( ...
            affineFamily{k} ...
            + affineFamily{k}.');
end

affineOptions = cfg.matrixFeatures;
affineOptions.name = 'affine family A_k - A_0';

affineReport = analyzeMatrixTrajectoryFeatures( ...
    affineFamily, affineOptions);

%% ========================================================================
% Positive-increment family: D_k^+
% ========================================================================

growthReport = analyzePositiveIncrementSubspace( ...
    A, x, cfg.growthSubspace);

positiveIncrementFamily = ...
    growthReport.positivePart;

positiveOptions = cfg.matrixFeatures;
positiveOptions.name = 'positive-increment family D_k^+';

positiveReport = analyzeMatrixTrajectoryFeatures( ...
    positiveIncrementFamily, positiveOptions);

%% ========================================================================
% Summary
% ========================================================================

reports = { ...
    linearReport, ...
    affineReport, ...
    positiveReport};

reportNames = { ...
    'Linear A_k', ...
    'Affine A_k-A_0', ...
    'Positive D_k^+'};

fprintf('\nCompact comparison:\n');

fprintf([ ...
    '%22s %8s %8s %8s %12s %12s\n'], ...
    'family', ...
    'rank1e-8', ...
    'rank1e-10', ...
    'rank1e-12', ...
    'm(spec<1e-2)', ...
    'm(spec<1e-3)');

for j = 1:numel(reports)
    report = reports{j};

    rank1e8 = getRankAtTolerance(report, 1e-8);
    rank1e10 = getRankAtTolerance(report, 1e-10);
    rank1e12 = getRankAtTolerance(report, 1e-12);

    dimension1e2 = getDimensionAtError(report, 1e-2);
    dimension1e3 = getDimensionAtError(report, 1e-3);

    fprintf([ ...
        '%22s %8d %8d %8d %12s %12s\n'], ...
        reportNames{j}, ...
        rank1e8, ...
        rank1e10, ...
        rank1e12, ...
        integerOrUnavailable(dimension1e2), ...
        integerOrUnavailable(dimension1e3));
end

%% ========================================================================
% Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Matrix feature spectra');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

colors = lines(3);

% -------------------------------------------------------------------------
% Singular values
% -------------------------------------------------------------------------

nexttile;
hold on;

for j = 1:3
    semilogy( ...
        1:numel(reports{j}.singularValues), ...
        max(reports{j}.singularValues, realmin), ...
        'LineWidth', 1.5, ...
        'Color', colors(j, :), ...
        'DisplayName', reportNames{j});
end

xlabel('Matrix feature index');
ylabel('Singular value');
title('Matrix-trajectory spectra');
legend('Location', 'best');
set(gca,'yscale','log')
grid on;
box on;

% -------------------------------------------------------------------------
% Cumulative Frobenius fraction
% -------------------------------------------------------------------------

nexttile;
hold on;

for j = 1:3
    plot( ...
        1:numel(reports{j}.cumulativeFrobeniusFraction), ...
        reports{j}.cumulativeFrobeniusFraction, ...
        'LineWidth', 1.5, ...
        'Color', colors(j, :), ...
        'DisplayName', reportNames{j});
end

xlabel('Matrix feature dimension');
ylabel('Cumulative Frobenius fraction');
title('Matrix-family compression');
legend('Location', 'best');
set(gca,'yscale','log')
grid on;
box on;
ylim([0, 1.001]);

% -------------------------------------------------------------------------
% Global spectral reconstruction error
% -------------------------------------------------------------------------

nexttile;
hold on;

for j = 1:3
    semilogy( ...
        reports{j}.candidateDimensions, ...
        max( ...
            reports{j}.maximumRelativeSpectralErrorGlobal, ...
            realmin), ...
        'o-', ...
        'LineWidth', 1.25, ...
        'Color', colors(j, :), ...
        'DisplayName', reportNames{j});
end

xlabel('Matrix feature dimension');
ylabel('Maximum spectral residual / global scale');
title('Global spectral reconstruction error');
legend('Location', 'best');
set(gca,'yscale','log')
grid on;
box on;

% -------------------------------------------------------------------------
% Local-normalized spectral reconstruction error
% -------------------------------------------------------------------------

nexttile;
hold on;

for j = 1:3
    semilogy( ...
        reports{j}.candidateDimensions, ...
        max( ...
            reports{j}.maximumRelativeSpectralErrorLocal, ...
            realmin), ...
        'o-', ...
        'LineWidth', 1.25, ...
        'Color', colors(j, :), ...
        'DisplayName', reportNames{j});
end

xlabel('Matrix feature dimension');
ylabel('Maximum local-relative spectral residual');
title('Worst local matrix reconstruction');
legend('Location', 'best');
set(gca,'yscale','log')
grid on;
box on;

sgtitle('Polyhedral-event matrix feature feasibility');

% -------------------------------------------------------------------------
% Coefficient trajectories for the first few affine features
% -------------------------------------------------------------------------

figure('Color', 'w', ...
    'Name', 'Affine matrix feature coefficients');

numCoefficientPlots = min( ...
    10, size(affineReport.coefficientTrajectories, 1));

hold on;

coefficientColors = lines(numCoefficientPlots);

for ell = 1:numCoefficientPlots
    coefficient = ...
        affineReport.coefficientTrajectories(ell, :);

    coefficientScale = max(abs(coefficient));

    if coefficientScale > 0
        coefficient = coefficient / coefficientScale;
    end

    plot( ...
        x, ...
        coefficient, ...
        'LineWidth', 1.25, ...
        'Color', coefficientColors(ell, :), ...
        'DisplayName', ...
        sprintf('feature %d', ell));
end

xlabel('Wall arc length');
ylabel('Normalized coefficient');
title('Leading affine matrix-feature trajectories');
legend('Location', 'eastoutside');
grid on;
box on;

%% ========================================================================
% Save
% ========================================================================

outputDirectory = fileparts(cfg.output.file);

if ~isempty(outputDirectory) && ...
        ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

save( ...
    cfg.output.file, ...
    'linearReport', ...
    'affineReport', ...
    'positiveReport', ...
    'growthReport', ...
    'cfg', ...
    'meta', ...
    '-v7.3');

fprintf('\nSaved matrix-feature report:\n  %s\n', ...
    cfg.output.file);

fprintf('\n============================================================\n');
fprintf('OWNS matrix-feature study completed.\n');
fprintf('============================================================\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function rankValue = getRankAtTolerance(report, tolerance)

[~, index] = min(abs( ...
    log10(report.rankTolerances) ...
    - log10(tolerance)));

rankValue = report.numericalRanks(index);

end

function dimension = getDimensionAtError(report, target)

[~, index] = min(abs( ...
    log10(report.spectralErrorTargets) ...
    - log10(target)));

dimension = ...
    report.dimensionForGlobalSpectralTarget(index);

end

function output = integerOrUnavailable(value)

if isnan(value)
    output = 'not reached';
else
    output = sprintf('%d', round(value));
end

end
