% runOWNSDeficitLoewnerEnvelopeStudy.m
%
% Compare two deterministic Loewner upper families:
%
%   1. Positive-variation envelope
%   2. Sequential deficit-based envelope
%
% No stochastic or angular sampling is used.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS deficit-based Loewner envelope study\n');
fprintf('============================================================\n\n');

%% ========================================================================
% Configuration
% ========================================================================

cfg = struct();

cfg.owns = struct();

cfg.owns.dataFile = ...
    ['/data2/jampauls/DataforFigures/SaveData/', ...
     'FirstPassingCode/data/ninth_Run_i_1_j_1_solution.mat'];

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
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.growthSubspace = struct();
cfg.growthSubspace.incrementScaling = 'perUnitX';

cfg.growthSubspace.positiveEigenvalueTolerances = ...
    [1e-6, 1e-8, 1e-10, 1e-12, 1e-14];

cfg.growthSubspace.positivePartTolerance = 1e-12;
cfg.growthSubspace.aggregateWeightMethod = 'magnitude';

cfg.growthSubspace.candidateRanks = ...
    [3, 5, 10, 15, 20, 30, 40, 60, 80];

cfg.growthSubspace.complementTolerance = 1e-10;
cfg.growthSubspace.verbose = false;

cfg.deficit = struct();
cfg.deficit.validationTolerance = 1e-10;
cfg.deficit.validateAllPrefixes = true;
cfg.deficit.verbose = true;

cfg.probability = struct();

cfg.probability.cfOptions = struct();
cfg.probability.cfOptions.absoluteTolerance = 1e-9;
cfg.probability.cfOptions.relativeTolerance = 1e-7;
cfg.probability.cfOptions.eigenvalueTolerance = 1e-13;
cfg.probability.cfOptions.equalEigenvalueTolerance = 1e-12;
cfg.probability.cfOptions.maxIntervalCount = 30000;
cfg.probability.cfOptions.warnOnLargeError = true;

cfg.output = struct();

cfg.output.file = fullfile( ...
    projectRoot, ...
    'results', ...
    'owns_deficit_loewner_envelope_study.mat');

cfg.output.saveMatrixFamilies = false;

%% ========================================================================
% Load matrix family
% ========================================================================

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x(:);
A = problem.A;
eThresh = problem.eThresh(:);

clear problem;

Nx = numel(x);
r = size(A{1}, 1);

thresholdScale = max(max(abs(eThresh)), 1);

if max(eThresh) - min(eThresh) > ...
        1e-12 * thresholdScale
    error( ...
        'runOWNSDeficitLoewnerEnvelopeStudy:NonconstantThreshold', ...
        'This study requires a constant threshold.');
end

threshold = eThresh(1);

fprintf('Problem dimensions:\n');
fprintf('  stations             = %d\n', Nx);
fprintf('  stochastic dimension = %d\n', r);
fprintf('  threshold            = %.12e\n\n', threshold);

%% ========================================================================
% Existing positive-variation family
% ========================================================================

fprintf('Building positive-variation family...\n');

growthReport = analyzePositiveIncrementSubspace( ...
    A, x, cfg.growthSubspace);

[positiveVariationFamily, positiveVariationDiagnostics] = ...
    buildPositiveVariationUpperFamily( ...
        A, x, growthReport);

%% ========================================================================
% Deficit-based family
% ========================================================================

fprintf('Building deficit-based family...\n');

[deficitFamily, deficitDiagnostics] = ...
    computeSequentialLoewnerEnvelope( ...
        A, cfg.deficit);

fprintf('\nDeficit-family matrix diagnostics:\n');
fprintf('  minimum dominance eigenvalue = %.12e\n', ...
    deficitDiagnostics.minimumDominanceEigenvalue);
fprintf('  minimum monotonicity eigenvalue = %.12e\n', ...
    deficitDiagnostics.minimumMonotonicityEigenvalue);
fprintf('  relative dominance defect = %.12e\n', ...
    deficitDiagnostics.relativeDominanceDefect);
fprintf('  relative monotonicity defect = %.12e\n', ...
    deficitDiagnostics.relativeMonotonicityDefect);
fprintf('  numerically dominating = %d\n', ...
    deficitDiagnostics.isNumericallyDominating);
fprintf('  numerically monotone   = %d\n\n', ...
    deficitDiagnostics.isNumericallyMonotone);

if ~deficitDiagnostics.isNumericallyDominating || ...
        ~deficitDiagnostics.isNumericallyMonotone
    error( ...
        'runOWNSDeficitLoewnerEnvelopeStudy:InvalidEnvelope', ...
        'The deficit-based family failed matrix-order validation.');
end

%% ========================================================================
% Deterministic scalar tails
% ========================================================================

localProbability = nan(Nx, 1);
positiveVariationProbability = nan(Nx, 1);
deficitProbability = nan(Nx, 1);

localErrorEstimate = zeros(Nx, 1);
positiveVariationErrorEstimate = zeros(Nx, 1);
deficitErrorEstimate = zeros(Nx, 1);

localWarning = cell(Nx, 1);
positiveVariationWarning = cell(Nx, 1);
deficitWarning = cell(Nx, 1);

fprintf('Evaluating deterministic quadratic-form tails...\n');

progressInterval = max(floor(Nx / 10), 1);

for k = 1:Nx
    Ak = real(0.5 * (A{k} + A{k}.'));

    Upv = real(0.5 * ( ...
        positiveVariationFamily{k} ...
        + positiveVariationFamily{k}.'));

    Ud = real(0.5 * ( ...
        deficitFamily{k} ...
        + deficitFamily{k}.'));

    lastwarn('');

    [localProbability(k), localInfo] = ...
        generalizedQuadraticFormTail( ...
            eig(Ak), threshold, 'real', ...
            cfg.probability.cfOptions);

    [localWarning{k}, ~] = lastwarn;

    localErrorEstimate(k) = ...
        extractIntegrationError(localInfo);

    lastwarn('');

    [positiveVariationProbability(k), pvInfo] = ...
        generalizedQuadraticFormTail( ...
            eig(Upv), threshold, 'real', ...
            cfg.probability.cfOptions);

    [positiveVariationWarning{k}, ~] = lastwarn;

    positiveVariationErrorEstimate(k) = ...
        extractIntegrationError(pvInfo);

    lastwarn('');

    [deficitProbability(k), deficitInfo] = ...
        generalizedQuadraticFormTail( ...
            eig(Ud), threshold, 'real', ...
            cfg.probability.cfOptions);

    [deficitWarning{k}, ~] = lastwarn;

    deficitErrorEstimate(k) = ...
        extractIntegrationError(deficitInfo);

    if mod(k, progressInterval) == 0 || k == Nx
        fprintf('  completed station %d of %d\n', k, Nx);
    end
end

%% ========================================================================
% Probability brackets and comparisons
% ========================================================================

lowerProbability = cummax(localProbability);

% Both matrix families are PSD-monotone, so their exact tails are
% monotone. Enforce only roundoff-level numerical monotonicity.
positiveVariationProbability = ...
    cummax(positiveVariationProbability);

deficitProbability = ...
    cummax(deficitProbability);

combinedUpperProbability = min( ...
    positiveVariationProbability, ...
    deficitProbability);

positiveVariationWidth = ...
    positiveVariationProbability - lowerProbability;

deficitWidth = ...
    deficitProbability - lowerProbability;

combinedWidth = ...
    combinedUpperProbability - lowerProbability;

probabilityImprovement = ...
    positiveVariationProbability ...
    - deficitProbability;

%% ========================================================================
% Matrix-size comparison
% ========================================================================

positiveVariationTrace = zeros(Nx, 1);
deficitTrace = zeros(Nx, 1);

positiveVariationSpectralNorm = zeros(Nx, 1);
deficitSpectralNorm = zeros(Nx, 1);

minimumPvMinusDeficitEigenvalue = zeros(Nx, 1);
minimumDeficitMinusPvEigenvalue = zeros(Nx, 1);

for k = 1:Nx
    Upv = positiveVariationFamily{k};
    Ud = deficitFamily{k};

    positiveVariationTrace(k) = trace(Upv);
    deficitTrace(k) = trace(Ud);

    positiveVariationSpectralNorm(k) = norm(Upv, 2);
    deficitSpectralNorm(k) = norm(Ud, 2);

    difference = Upv - Ud;
    difference = 0.5 * (difference + difference.');

    lambdaDifference = real(eig(difference));

    minimumPvMinusDeficitEigenvalue(k) = ...
        min(lambdaDifference);

    minimumDeficitMinusPvEigenvalue(k) = ...
        min(-lambdaDifference);
end

%% ========================================================================
% Summary
% ========================================================================

fprintf('\nDeterministic envelope comparison:\n');

fprintf('  terminal lower bound = %.12e\n', ...
    lowerProbability(end));

fprintf('  terminal positive-variation upper = %.12e\n', ...
    positiveVariationProbability(end));

fprintf('  terminal deficit upper            = %.12e\n', ...
    deficitProbability(end));

fprintf('  terminal combined upper           = %.12e\n', ...
    combinedUpperProbability(end));

fprintf('  terminal positive-variation width = %.12e\n', ...
    positiveVariationWidth(end));

fprintf('  terminal deficit width            = %.12e\n', ...
    deficitWidth(end));

fprintf('  terminal combined width           = %.12e\n', ...
    combinedWidth(end));

fprintf('  terminal deficit improvement      = %.12e\n', ...
    probabilityImprovement(end));

fprintf('\nEnvelope matrix comparison:\n');

fprintf('  terminal positive-variation trace = %.12e\n', ...
    positiveVariationTrace(end));

fprintf('  terminal deficit trace            = %.12e\n', ...
    deficitTrace(end));

fprintf('  terminal trace ratio deficit/PV   = %.12e\n', ...
    deficitTrace(end) ...
    / max(abs(positiveVariationTrace(end)), eps));

fprintf('  terminal positive-variation norm  = %.12e\n', ...
    positiveVariationSpectralNorm(end));

fprintf('  terminal deficit norm             = %.12e\n', ...
    deficitSpectralNorm(end));

fprintf('  min eig(PV-deficit) over x        = %.12e\n', ...
    min(minimumPvMinusDeficitEigenvalue));

fprintf('  min eig(deficit-PV) over x        = %.12e\n', ...
    min(minimumDeficitMinusPvEigenvalue));

fprintf('\nNumerical diagnostics:\n');

fprintf('  local-tail warnings = %d\n', ...
    nnz(~cellfun(@isempty, localWarning)));

fprintf('  positive-variation warnings = %d\n', ...
    nnz(~cellfun(@isempty, positiveVariationWarning)));

fprintf('  deficit-envelope warnings = %d\n', ...
    nnz(~cellfun(@isempty, deficitWarning)));

fprintf('  maximum local error estimate = %.12e\n', ...
    max(localErrorEstimate));

fprintf('  maximum PV error estimate    = %.12e\n', ...
    max(positiveVariationErrorEstimate));

fprintf('  maximum deficit error estimate = %.12e\n', ...
    max(deficitErrorEstimate));

fprintf('  minimum combined bracket margin = %.12e\n', ...
    min(combinedWidth));

%% ========================================================================
% Report
% ========================================================================

report = struct();

report.x = x;
report.threshold = threshold;

report.localProbability = localProbability;
report.lowerProbability = lowerProbability;

report.positiveVariationProbability = ...
    positiveVariationProbability;

report.deficitProbability = ...
    deficitProbability;

report.combinedUpperProbability = ...
    combinedUpperProbability;

report.positiveVariationWidth = ...
    positiveVariationWidth;

report.deficitWidth = deficitWidth;
report.combinedWidth = combinedWidth;

report.probabilityImprovement = ...
    probabilityImprovement;

report.positiveVariationTrace = ...
    positiveVariationTrace;

report.deficitTrace = deficitTrace;

report.positiveVariationSpectralNorm = ...
    positiveVariationSpectralNorm;

report.deficitSpectralNorm = ...
    deficitSpectralNorm;

report.minimumPvMinusDeficitEigenvalue = ...
    minimumPvMinusDeficitEigenvalue;

report.minimumDeficitMinusPvEigenvalue = ...
    minimumDeficitMinusPvEigenvalue;

report.localErrorEstimate = ...
    localErrorEstimate;

report.positiveVariationErrorEstimate = ...
    positiveVariationErrorEstimate;

report.deficitErrorEstimate = ...
    deficitErrorEstimate;

report.positiveVariationDiagnostics = ...
    positiveVariationDiagnostics;

report.deficitDiagnostics = ...
    deficitDiagnostics;

report.localWarning = localWarning;
report.positiveVariationWarning = ...
    positiveVariationWarning;

report.deficitWarning = deficitWarning;

report.isDeterministic = true;
report.usesRandomSampling = false;

report.interpretation = ...
    ['Comparison of positive-variation and sequential deficit-based ', ...
     'Loewner envelopes. The smaller valid scalar tail provides the ', ...
     'combined deterministic upper bound.'];

%% ========================================================================
% Plots
% ========================================================================

figure( ...
    'Color', 'w', ...
    'Name', 'Deficit-based Loewner envelope comparison');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
hold on;

plot(x, lowerProbability, ...
    'k--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Local lower bound');

plot(x, positiveVariationProbability, ...
    ':', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Positive-variation upper');

plot(x, deficitProbability, ...
    '-', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Deficit upper');

plot(x, combinedUpperProbability, ...
    '-.', ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Combined upper');

xlabel('Wall arc length');
ylabel('Probability');
title('Deterministic probability bounds');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

plot(x, positiveVariationWidth, ...
    ':', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Positive variation');

plot(x, deficitWidth, ...
    '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Deficit');

plot(x, combinedWidth, ...
    '-.', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Combined');

xlabel('Wall arc length');
ylabel('Upper bound - lower bound');
title('Bracket widths');
legend('Location', 'best');
grid on;
box on;

nexttile;
hold on;

plot(x, positiveVariationTrace, ...
    ':', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Positive variation');

plot(x, deficitTrace, ...
    '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Deficit');

xlabel('Wall arc length');
ylabel('Envelope trace');
title('Envelope size');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, probabilityImprovement, ...
    'LineWidth', 1.5);

yline(0, 'k--');

xlabel('Wall arc length');
ylabel('PV upper - deficit upper');
title('Deficit-envelope probability improvement');
grid on;
box on;

sgtitle('Deterministic Loewner-envelope comparison');

%% ========================================================================
% Save
% ========================================================================

outputDirectory = fileparts(cfg.output.file);

if ~isempty(outputDirectory) && ...
        ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

if cfg.output.saveMatrixFamilies
    save( ...
        cfg.output.file, ...
        'report', ...
        'cfg', ...
        'positiveVariationFamily', ...
        'deficitFamily', ...
        'A', ...
        '-v7.3');
else
    save( ...
        cfg.output.file, ...
        'report', ...
        'cfg', ...
        '-v7.3');
end

fprintf('\nSaved deficit-envelope study:\n');
fprintf('  %s\n', cfg.output.file);

fprintf('\n============================================================\n');
fprintf('OWNS deficit-based Loewner envelope study completed.\n');
fprintf('============================================================\n');

% =========================================================================
function errorEstimate = extractIntegrationError(info)

errorEstimate = 0;

if isstruct(info) && ...
        isfield(info, 'integrationErrorEstimate') && ...
        ~isempty(info.integrationErrorEstimate) && ...
        isfinite(info.integrationErrorEstimate)

    errorEstimate = max( ...
        real(info.integrationErrorEstimate), 0);
end

end
