% runOWNSMatrixFeatureProbabilityStudy.m
%
% Validate matrix-feature truncations against the full first-transition CDF
% using common scrambled-Sobol RQMC directions.
%
% RQMC is used only as a validation instrument in this study.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS matrix-feature first-passage validation\n');
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

% Use the same provisional threshold as the previous probability studies.
cfg.threshold.method = 'auto';
cfg.threshold.autoMaxMeanFactor = 1;
cfg.threshold.warnIfProvisional = false;

cfg.numerics = struct();
cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = false;
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
cfg.angular.randomSeed = 246802;
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;
cfg.angular.scramble = true;
cfg.angular.probabilityClip = 1e-12;

cfg.featureStudy = struct();

cfg.featureStudy.dimensions = ...
    [3, 5, 10, 15, 20, 30];

cfg.featureStudy.savedFeatureReport = fullfile( ...
    projectRoot, ...
    'results', ...
    'owns_polyhedral_structure_study.mat');

cfg.featureStudy.numBasisMatricesToStore = ...
    max(cfg.featureStudy.dimensions);

cfg.output = struct();

cfg.output.file = fullfile( ...
    projectRoot, ...
    'results', ...
    'owns_matrix_feature_probability_study.mat');

%% ========================================================================
% Load OWNS matrix family
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x;
A = problem.A;
eThresh = problem.eThresh;
meta = problem.meta;

clear problem;

if max(eThresh) - min(eThresh) > ...
        1e-12 * max(abs(eThresh))
    error('runOWNSMatrixFeatureProbabilityStudy:Threshold', ...
        'This study currently requires a constant threshold.');
end

threshold = eThresh(1);

%% ========================================================================
% Load or build direct-SVD feature report
% ========================================================================

featureReport = [];

if isfile(cfg.featureStudy.savedFeatureReport)
    saved = load( ...
        cfg.featureStudy.savedFeatureReport, ...
        'featureReport');

    if isfield(saved, 'featureReport')
        featureReport = saved.featureReport;

        fprintf('Loaded existing direct-SVD feature report:\n');
        fprintf('  %s\n\n', ...
            cfg.featureStudy.savedFeatureReport);
    end
end

requiredFeatureDimension = ...
    max(cfg.featureStudy.dimensions);

needsRebuild = isempty(featureReport) ...
    || ~isfield(featureReport, 'basisMatrices') ...
    || featureReport.numBasisMatricesStored ...
        < requiredFeatureDimension ...
    || ~isfield(featureReport, 'coefficientTrajectories');

if needsRebuild
    fprintf('Recomputing direct matrix-trajectory SVD...\n');

    svdOptions = struct();

    svdOptions.name = 'linear family A_k';

    svdOptions.candidateDimensions = [ ...
        0, ...
        cfg.featureStudy.dimensions, ...
        40, 60, 80, 100];

    svdOptions.rankTolerances = ...
        [1e-6, 1e-8, 1e-10, 1e-12, 1e-14];

    svdOptions.spectralErrorTargets = ...
        [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4, 1e-6];

    svdOptions.numBasisMatricesToStore = ...
        cfg.featureStudy.numBasisMatricesToStore;

    svdOptions.verbose = true;

    featureReport = ...
        analyzeMatrixTrajectoryFeaturesDirectSVD( ...
            A, svdOptions);
end

%% ========================================================================
% Coupled probability study
% ========================================================================

study = computeCoupledMatrixFeatureProbabilityStudy( ...
    x, ...
    A, ...
    featureReport, ...
    cfg.featureStudy.dimensions, ...
    threshold, ...
    cfg.angular, ...
    cfg.maxDetection, ...
    cfg.numerics);

%% ========================================================================
% Print summary
% ========================================================================

fprintf('\nMatrix-feature first-passage summary:\n');

fprintf([ ...
    '%6s %14s %12s %14s %14s %14s ', ...
    '%12s %12s\n'], ...
    'm', ...
    'max CDF diff', ...
    'paired SE', ...
    'terminal diff', ...
    'max local diff', ...
    'false positive', ...
    'false neg.', ...
    'PSD defect');

for iFeature = 1:numel(study.featureDimensions)
    diagnostics = ...
        study.matrixDiagnostics{iFeature};

    fprintf([ ...
        '%6d %14.6e %12.4e %14.6e %14.6e ', ...
        '%14.6e %12.6e %12.4e\n'], ...
        study.featureDimensions(iFeature), ...
        study.maximumAbsoluteCdfDifference(iFeature), ...
        study.maximumAbsoluteCdfDifferenceStandardError(iFeature), ...
        study.terminalCdfDifference(iFeature), ...
        study.maximumAbsoluteLocalDifference(iFeature), ...
        study.maximumFalsePositiveProbability(iFeature), ...
        study.maximumFalseNegativeProbability(iFeature), ...
        diagnostics.maximumRelativePsdDefect);
end

fprintf('\nFull-model terminal probability:\n');
fprintf('  F(xmax) = %.12e\n', study.Ffull(end));
fprintf('  SE      = %.12e\n\n', ...
    study.FfullStandardError(end));

%% ========================================================================
% Plots
% ========================================================================

featureDimensions = study.featureDimensions;
numFeatures = numel(featureDimensions);

colors = lines(numFeatures);

figure('Color', 'w', ...
    'Name', 'Matrix-feature first-passage accuracy');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% Maximum CDF error
% -------------------------------------------------------------------------

nexttile;

lowerError = min( ...
    2 * study.maximumAbsoluteCdfDifferenceStandardError, ...
    0.95 * study.maximumAbsoluteCdfDifference);

upperError = ...
    2 * study.maximumAbsoluteCdfDifferenceStandardError;

errorbar( ...
    featureDimensions, ...
    study.maximumAbsoluteCdfDifference, ...
    lowerError, ...
    upperError, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

set(gca, 'YScale', 'log');

xlabel('Matrix feature dimension m');
ylabel('Maximum absolute CDF difference');
title('First-passage feature-truncation error');
grid on;
box on;

% -------------------------------------------------------------------------
% False-positive and false-negative probabilities
% -------------------------------------------------------------------------

nexttile;
hold on;

semilogy( ...
    featureDimensions, ...
    max(study.maximumFalsePositiveProbability, realmin), ...
    'o-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'False positive');

semilogy( ...
    featureDimensions, ...
    max(study.maximumFalseNegativeProbability, realmin), ...
    's--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'False negative');

xlabel('Matrix feature dimension m');
ylabel('Maximum disagreement probability');
title('Transition classification disagreement');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Complete CDFs
% -------------------------------------------------------------------------

nexttile;
hold on;

plot(x, study.Ffull, ...
    'k', ...
    'LineWidth', 2, ...
    'DisplayName', 'Full matrix family');

for iFeature = 1:numFeatures
    plot( ...
        x, ...
        study.Ffeature(:, iFeature), ...
        'LineWidth', 1.15, ...
        'Color', colors(iFeature, :), ...
        'DisplayName', ...
        sprintf('m=%d', ...
            featureDimensions(iFeature)));
end

xlabel('Wall arc length');
ylabel('F_{X_{tr}}(x)');
title('Full and feature-truncated CDFs');
legend('Location', 'eastoutside');
grid on;
box on;

% -------------------------------------------------------------------------
% CDF difference over x
% -------------------------------------------------------------------------

nexttile;
hold on;

for iFeature = 1:numFeatures
    plot( ...
        x, ...
        study.cdfDifference(:, iFeature), ...
        'LineWidth', 1.2, ...
        'Color', colors(iFeature, :), ...
        'DisplayName', ...
        sprintf('m=%d', ...
            featureDimensions(iFeature)));
end

xlabel('Wall arc length');
ylabel('F_m-F_{full}');
title('Paired first-passage error');
legend('Location', 'best');
grid on;
box on;

sgtitle('Matrix-feature probability validation');

figure('Color', 'w', ...
    'Name', 'Matrix-feature local and envelope errors');

tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
hold on;

for iFeature = 1:numFeatures
    plot( ...
        x, ...
        study.localDifference(:, iFeature), ...
        'LineWidth', 1.2, ...
        'Color', colors(iFeature, :), ...
        'DisplayName', ...
        sprintf('m=%d', ...
            featureDimensions(iFeature)));
end

xlabel('Wall arc length');
ylabel('p_{local,m}-p_{local,full}');
title('Local probability error');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

for iFeature = 1:numFeatures
    semilogy( ...
        x, ...
        max( ...
            study.weightedEnvelopeAbsoluteError(:, iFeature), ...
            realmin), ...
        'LineWidth', 1.2, ...
        'Color', colors(iFeature, :), ...
        'DisplayName', ...
        sprintf('m=%d', ...
            featureDimensions(iFeature)));
end

xlabel('Wall arc length');
ylabel('Weighted mean envelope error');
title('Directional historical-envelope error');
legend('Location', 'best');
grid on;
box on;

sgtitle('Matrix-feature process diagnostics');

%% ========================================================================
% Save
% ========================================================================

outputDirectory = fileparts(cfg.output.file);

if ~isempty(outputDirectory) && ...
        ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

% The reconstructed matrix families can make the saved structure large.
% Remove them before saving unless they are specifically needed later.
studyToSave = study;
studyToSave.Afeature = [];

save( ...
    cfg.output.file, ...
    'studyToSave', ...
    'featureReport', ...
    'cfg', ...
    'meta', ...
    '-v7.3');

fprintf('Saved matrix-feature probability study:\n');
fprintf('  %s\n\n', cfg.output.file);

fprintf('============================================================\n');
fprintf('OWNS matrix-feature probability study completed.\n');
fprintf('============================================================\n');
