% runOWNSDeterministicLoewnerBoundStudy.m
%
% Deterministically evaluate first-transition lower and upper bounds:
%
%   max_{j<=k} P(w.'*A_j*w >= c)
%       <= P(max_{j<=k} w.'*A_j*w >= c)
%       <= P(w.'*U_k*w >= c),
%
% where U_k is the positive-variation Loewner upper family.
%
% No Gaussian vectors, angular directions, Sobol' points, or other random
% samples are used.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS deterministic Loewner probability-bound study\n');
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

% These ranks are not used by the deterministic envelope calculation, but
% analyzePositiveIncrementSubspace may expect this configuration field.
cfg.growthSubspace.candidateRanks = ...
    [3, 5, 10, 15, 20, 30, 40, 60, 80];

cfg.growthSubspace.complementTolerance = 1e-10;
cfg.growthSubspace.verbose = false;

cfg.probability = struct();

cfg.probability.cfOptions = struct();
cfg.probability.cfOptions.absoluteTolerance = 1e-9;
cfg.probability.cfOptions.relativeTolerance = 1e-7;
cfg.probability.cfOptions.eigenvalueTolerance = 1e-13;
cfg.probability.cfOptions.equalEigenvalueTolerance = 1e-12;
cfg.probability.cfOptions.maxIntervalCount = 30000;

% Keep the project's own diagnostic enabled. Unexpected quadgk warnings
% will also remain visible rather than being suppressed.
cfg.probability.cfOptions.warnOnLargeError = true;

% Estimated CF integration errors are used to form conservative numerical
% curves. They are not mathematically rigorous interval enclosures.
cfg.probability.useEstimatedErrorMargins = true;

cfg.output = struct();

cfg.output.file = fullfile( ...
    projectRoot, ...
    'results', ...
    'owns_deterministic_loewner_bound_study.mat');

% Storing the complete full-size matrix families substantially increases
% the output file size. The default report stores only their spectra and
% diagnostics.
cfg.output.saveMatrixFamilies = false;

%% ========================================================================
% Load OWNS matrix family
% ========================================================================

fprintf('Loading OWNS problem:\n');
fprintf('  %s\n\n', cfg.owns.dataFile);

problem = loadOWNSProblem( ...
    cfg.owns.dataFile, cfg);

x = problem.x(:);
A = problem.A;
eThresh = problem.eThresh(:);

clear problem;

Nx = numel(x);
r = size(A{1}, 1);

if numel(A) ~= Nx || numel(eThresh) ~= Nx
    error( ...
        'runOWNSDeterministicLoewnerBoundStudy:LengthMismatch', ...
        'The coordinate, matrix, and threshold families must agree.');
end

thresholdScale = max(max(abs(eThresh)), 1);

if max(eThresh) - min(eThresh) > ...
        1e-12 * thresholdScale
    error( ...
        'runOWNSDeterministicLoewnerBoundStudy:NonconstantThreshold', ...
        'This study currently requires a constant threshold.');
end

threshold = eThresh(1);

fprintf('Problem dimensions:\n');
fprintf('  stations             = %d\n', Nx);
fprintf('  stochastic dimension = %d\n', r);
fprintf('  threshold            = %.12e\n\n', threshold);

%% ========================================================================
% Build the existing positive-variation upper family
% ========================================================================

fprintf('Analyzing positive increments...\n');

growthReport = analyzePositiveIncrementSubspace( ...
    A, x, cfg.growthSubspace);

fprintf('Building positive-variation upper family...\n');

[Ufamily, upperDiagnostics] = ...
    buildPositiveVariationUpperFamily( ...
        A, x, growthReport);

fprintf('\nMatrix-order diagnostics:\n');

fprintf('  minimum dominance eigenvalue = %.12e\n', ...
    upperDiagnostics.minimumDominanceEigenvalue);

fprintf('  minimum monotone increment eigenvalue = %.12e\n', ...
    upperDiagnostics.minimumMonotonicIncrementEigenvalue);

fprintf('  relative dominance defect = %.12e\n', ...
    upperDiagnostics.relativeDominanceDefect);

fprintf('  relative monotonicity defect = %.12e\n', ...
    upperDiagnostics.relativeMonotonicityDefect);

fprintf('  numerically dominating = %d\n', ...
    upperDiagnostics.isNumericallyDominating);

fprintf('  numerically monotone   = %d\n\n', ...
    upperDiagnostics.isNumericallyMonotone);

if ~upperDiagnostics.isNumericallyDominating
    error( ...
        'runOWNSDeterministicLoewnerBoundStudy:InvalidEnvelope', ...
        ['The positive-variation family did not pass the direct ', ...
         'Loewner-dominance check.']);
end

%% ========================================================================
% Deterministic local and envelope probabilities
% ========================================================================

fprintf('Evaluating deterministic scalar quadratic-form tails...\n');

localProbability = nan(Nx, 1);
upperProbability = nan(Nx, 1);

localErrorEstimate = zeros(Nx, 1);
upperErrorEstimate = zeros(Nx, 1);

localMethod = cell(Nx, 1);
upperMethod = cell(Nx, 1);

localMinimumEigenvalue = nan(Nx, 1);
upperMinimumEigenvalue = nan(Nx, 1);

localMaximumEigenvalue = nan(Nx, 1);
upperMaximumEigenvalue = nan(Nx, 1);

localTrace = nan(Nx, 1);
upperTrace = nan(Nx, 1);

localWarningMessage = cell(Nx, 1);
localWarningIdentifier = cell(Nx, 1);

upperWarningMessage = cell(Nx, 1);
upperWarningIdentifier = cell(Nx, 1);

progressInterval = max(floor(Nx / 10), 1);

for k = 1:Nx
    Ak = 0.5 * (A{k} + A{k}.');
    Uk = 0.5 * (Ufamily{k} + Ufamily{k}.');

    lambdaA = real(eig(Ak));
    lambdaU = real(eig(Uk));

    localMinimumEigenvalue(k) = min(lambdaA);
    upperMinimumEigenvalue(k) = min(lambdaU);

    localMaximumEigenvalue(k) = max(lambdaA);
    upperMaximumEigenvalue(k) = max(lambdaU);

    localTrace(k) = sum(lambdaA);
    upperTrace(k) = sum(lambdaU);

    % ---------------------------------------------------------------------
    % Local event probability
    % ---------------------------------------------------------------------

    lastwarn('');

    [localProbability(k), localInfo] = ...
        generalizedQuadraticFormTail( ...
            lambdaA, ...
            threshold, ...
            'real', ...
            cfg.probability.cfOptions);

    [localWarningMessage{k}, ...
        localWarningIdentifier{k}] = lastwarn;

    localErrorEstimate(k) = ...
        extractIntegrationError(localInfo);

    localMethod{k} = ...
        extractMethod(localInfo);

    % ---------------------------------------------------------------------
    % Envelope event probability
    % ---------------------------------------------------------------------

    lastwarn('');

    [upperProbability(k), upperInfo] = ...
        generalizedQuadraticFormTail( ...
            lambdaU, ...
            threshold, ...
            'real', ...
            cfg.probability.cfOptions);

    [upperWarningMessage{k}, ...
        upperWarningIdentifier{k}] = lastwarn;

    upperErrorEstimate(k) = ...
        extractIntegrationError(upperInfo);

    upperMethod{k} = ...
        extractMethod(upperInfo);

    if mod(k, progressInterval) == 0 || k == Nx
        fprintf('  completed station %d of %d\n', k, Nx);
    end
end

%% ========================================================================
% Form deterministic first-transition brackets
% ========================================================================

% Central numerical estimates.
lowerProbability = cummax(localProbability);

% Because U_k is Loewner-monotone, its exact tail should be monotone.
% cummax prevents tiny inversion fluctuations from violating this known
% analytical property.
upperProbabilityMonotone = cummax(upperProbability);

if cfg.probability.useEstimatedErrorMargins
    localProbabilityConservative = max( ...
        localProbability - localErrorEstimate, 0);

    envelopeProbabilityConservative = min( ...
        upperProbability + upperErrorEstimate, 1);
else
    localProbabilityConservative = localProbability;
    envelopeProbabilityConservative = upperProbability;
end

lowerProbabilityConservative = ...
    cummax(localProbabilityConservative);

upperProbabilityConservative = ...
    cummax(envelopeProbabilityConservative);

centralBracketWidth = ...
    upperProbabilityMonotone - lowerProbability;

conservativeBracketWidth = ...
    upperProbabilityConservative ...
    - lowerProbabilityConservative;

centralValidityMargin = ...
    upperProbabilityMonotone - lowerProbability;

conservativeValidityMargin = ...
    upperProbabilityConservative ...
    - lowerProbabilityConservative;

%% ========================================================================
% Numerical diagnostics
% ========================================================================

numLocalWarnings = nnz( ...
    ~cellfun(@isempty, localWarningMessage));

numUpperWarnings = nnz( ...
    ~cellfun(@isempty, upperWarningMessage));

maximumLocalErrorEstimate = ...
    max(localErrorEstimate);

maximumUpperErrorEstimate = ...
    max(upperErrorEstimate);

maximumUpperMonotonicityDrop = max([ ...
    0; ...
    upperProbability(1:end-1) ...
        - upperProbability(2:end)]);

minimumCentralBracketMargin = ...
    min(centralValidityMargin);

minimumConservativeBracketMargin = ...
    min(conservativeValidityMargin);

if minimumCentralBracketMargin < -1e-6
    warning( ...
        'runOWNSDeterministicLoewnerBoundStudy:InvalidCentralBracket', ...
        ['The central numerical upper curve falls below the local ', ...
         'lower curve by %.6e. Inspect CF-inversion accuracy.'], ...
        -minimumCentralBracketMargin);
end

%% ========================================================================
% Assemble report
% ========================================================================

report = struct();

report.x = x;
report.threshold = threshold;

report.localProbability = localProbability;
report.lowerProbability = lowerProbability;

report.upperProbability = upperProbability;
report.upperProbabilityMonotone = ...
    upperProbabilityMonotone;

report.localErrorEstimate = localErrorEstimate;
report.upperErrorEstimate = upperErrorEstimate;

report.lowerProbabilityConservative = ...
    lowerProbabilityConservative;

report.upperProbabilityConservative = ...
    upperProbabilityConservative;

report.centralBracketWidth = ...
    centralBracketWidth;

report.conservativeBracketWidth = ...
    conservativeBracketWidth;

report.localMethod = localMethod;
report.upperMethod = upperMethod;

report.localMinimumEigenvalue = ...
    localMinimumEigenvalue;

report.upperMinimumEigenvalue = ...
    upperMinimumEigenvalue;

report.localMaximumEigenvalue = ...
    localMaximumEigenvalue;

report.upperMaximumEigenvalue = ...
    upperMaximumEigenvalue;

report.localTrace = localTrace;
report.upperTrace = upperTrace;

report.localWarningMessage = ...
    localWarningMessage;

report.localWarningIdentifier = ...
    localWarningIdentifier;

report.upperWarningMessage = ...
    upperWarningMessage;

report.upperWarningIdentifier = ...
    upperWarningIdentifier;

report.numLocalWarnings = numLocalWarnings;
report.numUpperWarnings = numUpperWarnings;

report.maximumLocalErrorEstimate = ...
    maximumLocalErrorEstimate;

report.maximumUpperErrorEstimate = ...
    maximumUpperErrorEstimate;

report.maximumUpperMonotonicityDrop = ...
    maximumUpperMonotonicityDrop;

report.minimumCentralBracketMargin = ...
    minimumCentralBracketMargin;

report.minimumConservativeBracketMargin = ...
    minimumConservativeBracketMargin;

report.upperDiagnostics = upperDiagnostics;

report.isDeterministic = true;
report.usesRandomSampling = false;

report.interpretation = [ ...
    'The lower curve is the running maximum of deterministic local ', ...
    'quadratic-form probabilities. The upper curve is the deterministic ', ...
    'quadratic-form probability of the positive-variation Loewner ', ...
    'envelope. Numerical CF integration errors are estimates rather ', ...
    'than formally verified interval bounds.'];

%% ========================================================================
% Print summary
% ========================================================================

[maximumLocalProbability, maximumLocalIndex] = ...
    max(localProbability);

[maximumCentralWidth, maximumCentralWidthIndex] = ...
    max(centralBracketWidth);

[maximumConservativeWidth, ...
    maximumConservativeWidthIndex] = ...
    max(conservativeBracketWidth);

fprintf('\nDeterministic probability-bound summary:\n');

fprintf('  maximum local probability = %.12e at station %d\n', ...
    maximumLocalProbability, maximumLocalIndex);

fprintf('  terminal lower bound      = %.12e\n', ...
    lowerProbability(end));

fprintf('  terminal upper bound      = %.12e\n', ...
    upperProbabilityMonotone(end));

fprintf('  terminal bracket width    = %.12e\n', ...
    centralBracketWidth(end));

fprintf('  conservative terminal lower = %.12e\n', ...
    lowerProbabilityConservative(end));

fprintf('  conservative terminal upper = %.12e\n', ...
    upperProbabilityConservative(end));

fprintf('  conservative terminal width = %.12e\n', ...
    conservativeBracketWidth(end));

fprintf('\nComplete-curve diagnostics:\n');

fprintf('  maximum central bracket width = %.12e at station %d\n', ...
    maximumCentralWidth, maximumCentralWidthIndex);

fprintf('  maximum conservative width    = %.12e at station %d\n', ...
    maximumConservativeWidth, maximumConservativeWidthIndex);

fprintf('  minimum central bracket margin = %.12e\n', ...
    minimumCentralBracketMargin);

fprintf('  minimum conservative margin    = %.12e\n', ...
    minimumConservativeBracketMargin);

fprintf('  maximum upper monotonicity drop = %.12e\n', ...
    maximumUpperMonotonicityDrop);

fprintf('\nScalar inversion diagnostics:\n');

fprintf('  local-tail warnings    = %d of %d\n', ...
    numLocalWarnings, Nx);

fprintf('  envelope-tail warnings = %d of %d\n', ...
    numUpperWarnings, Nx);

fprintf('  maximum local error estimate = %.12e\n', ...
    maximumLocalErrorEstimate);

fprintf('  maximum upper error estimate = %.12e\n', ...
    maximumUpperErrorEstimate);

fprintf('\nEnvelope growth diagnostics:\n');

fprintf('  initial envelope trace = %.12e\n', ...
    upperTrace(1));

fprintf('  terminal envelope trace = %.12e\n', ...
    upperTrace(end));

fprintf('  trace growth factor = %.12e\n', ...
    upperTrace(end) / max(abs(upperTrace(1)), eps));

fprintf('  terminal envelope spectral norm = %.12e\n', ...
    upperMaximumEigenvalue(end));

%% ========================================================================
% Plot
% ========================================================================

figure( ...
    'Color', 'w', ...
    'Name', 'Deterministic Loewner probability bounds');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% Deterministic probability bracket
% -------------------------------------------------------------------------

nexttile;
hold on;

plot( ...
    x, ...
    lowerProbability, ...
    '--', ...
    'LineWidth', 1.6, ...
    'DisplayName', 'Local-event lower bound');

plot( ...
    x, ...
    upperProbabilityMonotone, ...
    ':', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Loewner-envelope upper bound');

xlabel('Wall arc length');
ylabel('Probability');
title('Deterministic first-transition bracket');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Bracket width
% -------------------------------------------------------------------------

nexttile;

plot( ...
    x, ...
    centralBracketWidth, ...
    'LineWidth', 1.5);

xlabel('Wall arc length');
ylabel('Upper bound - lower bound');
title('Deterministic bracket width');
grid on;
box on;

% -------------------------------------------------------------------------
% Matrix traces
% -------------------------------------------------------------------------

nexttile;
hold on;

plot( ...
    x, ...
    localTrace, ...
    'LineWidth', 1.4, ...
    'DisplayName', 'trace(A_k)');

plot( ...
    x, ...
    upperTrace, ...
    'LineWidth', 1.4, ...
    'DisplayName', 'trace(U_k)');

xlabel('Wall arc length');
ylabel('Matrix trace');
title('Positive-variation accumulation');
legend('Location', 'best');
grid on;
box on;

% -------------------------------------------------------------------------
% Numerical inversion error estimates
% -------------------------------------------------------------------------

nexttile;
hold on;

semilogy( ...
    x, ...
    max(localErrorEstimate, realmin), ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Local tails');

semilogy( ...
    x, ...
    max(upperErrorEstimate, realmin), ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Envelope tails');

xlabel('Wall arc length');
ylabel('Estimated probability error');
title('Scalar inversion diagnostics');
legend('Location', 'best');
grid on;
box on;

sgtitle('Deterministic positive-variation Loewner bounds');

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
        'Ufamily', ...
        'A', ...
        '-v7.3');
else
    save( ...
        cfg.output.file, ...
        'report', ...
        'cfg', ...
        '-v7.3');
end

fprintf('\nSaved deterministic Loewner-bound study:\n');
fprintf('  %s\n', cfg.output.file);

fprintf('\n============================================================\n');
fprintf('OWNS deterministic Loewner-bound study completed.\n');
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

% =========================================================================
function method = extractMethod(info)

method = 'unknown';

if isstruct(info) && ...
        isfield(info, 'method') && ...
        ~isempty(info.method)

    method = char(info.method);
end

end