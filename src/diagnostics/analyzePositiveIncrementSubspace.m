function report = analyzePositiveIncrementSubspace(A, x, options)
%ANALYZEPOSITIVEINCREMENTSUBSPACE Analyze directions of energy growth.
%
% For each interval:
%
%   D_k = (A_{k+1}-A_k)/(x_{k+1}-x_k)
%
% and:
%
%   D_k = D_k^+ - D_k^-,
%
% where D_k^+ and D_k^- are positive semidefinite.
%
% A fixed growth basis is constructed from an aggregate of D_k^+.
% Candidate basis dimensions are tested for:
%
%   positive-increment reconstruction;
%   trace capture;
%   monotonicity of the orthogonal complement;
%   active-complement coupling.
%
% IMPORTANT:
%   These are floating-point diagnostics. Numerical semidefiniteness is
%   not a symbolic exact certificate.
%
% Inputs:
%   A        cell array of real symmetric r-by-r matrices
%   x        Nx-by-1 strictly increasing coordinate
%   options  diagnostic configuration
%
% MATLAB version: R2020b

if nargin < 3
    options = struct();
end

x = x(:);
Nx = numel(x);

if numel(A) ~= Nx
    error('analyzePositiveIncrementSubspace:LengthMismatch', ...
        'A must contain one matrix per streamwise station.');
end

if Nx < 2
    error('analyzePositiveIncrementSubspace:InsufficientStations', ...
        'At least two streamwise stations are required.');
end

if any(~isfinite(x)) || any(diff(x) <= 0)
    error('analyzePositiveIncrementSubspace:InvalidCoordinate', ...
        'x must be finite and strictly increasing.');
end

r = size(A{1}, 1);
numIntervals = Nx - 1;

incrementScaling = lower(char(getOption( ...
    options, 'incrementScaling', 'perUnitX')));

positiveEigenvalueTolerances = getOption( ...
    options, ...
    'positiveEigenvalueTolerances', ...
    [1e-6, 1e-8, 1e-10, 1e-12, 1e-14]);

positiveEigenvalueTolerances = ...
    positiveEigenvalueTolerances(:).';

positivePartTolerance = getOption( ...
    options, 'positivePartTolerance', 1e-12);

aggregateWeightMethod = lower(char(getOption( ...
    options, 'aggregateWeightMethod', 'magnitude')));

candidateRanks = getOption( ...
    options, ...
    'candidateRanks', ...
    [1, 2, 3, 5, 10, 20, 30, 40, 60, 80, 100]);

candidateRanks = unique(round(candidateRanks(:).'));
candidateRanks = candidateRanks( ...
    candidateRanks >= 1 & candidateRanks <= r);

complementTolerance = getOption( ...
    options, 'complementTolerance', 1e-10);

verbose = getOption(options, 'verbose', true);

dx = diff(x);

% -------------------------------------------------------------------------
% Construct all matrix increments before selecting a positive threshold.
% -------------------------------------------------------------------------

increment = cell(numIntervals, 1);

incrementEigenvalues = zeros(r, numIntervals);
incrementSpectralNorm = zeros(numIntervals, 1);
incrementFrobeniusNorm = zeros(numIntervals, 1);

for k = 1:numIntervals
    Ak = 0.5 * (A{k} + A{k}.');
    Akp1 = 0.5 * (A{k + 1} + A{k + 1}.');

    if any(size(Ak) ~= [r, r]) || ...
            any(size(Akp1) ~= [r, r])
        error('analyzePositiveIncrementSubspace:DimensionMismatch', ...
            'The A matrix dimensions are inconsistent.');
    end

    Dk = Akp1 - Ak;

    switch incrementScaling
        case 'perunitx'
            Dk = Dk / dx(k);

        case 'raw'
            % Use A_{k+1}-A_k without spatial normalization.

        otherwise
            error('analyzePositiveIncrementSubspace:UnknownScaling', ...
                'Unknown incrementScaling "%s".', ...
                options.incrementScaling);
    end

    Dk = 0.5 * (Dk + Dk.');

    increment{k} = Dk;

    lambda = sort(real(eig(Dk)), 'descend');

    incrementEigenvalues(:, k) = lambda;
    incrementSpectralNorm(k) = max(abs(lambda));
    incrementFrobeniusNorm(k) = norm(Dk, 'fro');
end

referenceIncrementScale = max(incrementSpectralNorm);

if referenceIncrementScale <= 0
    referenceIncrementScale = 1;
end

% -------------------------------------------------------------------------
% Split every increment into positive and negative PSD parts.
% -------------------------------------------------------------------------

positivePart = cell(numIntervals, 1);
negativePart = cell(numIntervals, 1);

positiveRank = zeros( ...
    numIntervals, numel(positiveEigenvalueTolerances));

positiveTrace = zeros(numIntervals, 1);
negativeTrace = zeros(numIntervals, 1);

positiveSpectralNorm = zeros(numIntervals, 1);
negativeSpectralNorm = zeros(numIntervals, 1);

positiveFrobeniusNorm = zeros(numIntervals, 1);
negativeFrobeniusNorm = zeros(numIntervals, 1);

retainedPositiveRank = zeros(numIntervals, 1);

for k = 1:numIntervals
    Dk = increment{k};

    [V, D] = eig(Dk);
    lambda = real(diag(D));

    [lambda, order] = sort(lambda, 'descend');
    V = real(V(:, order));

    for iTol = 1:numel(positiveEigenvalueTolerances)
        threshold = ...
            positiveEigenvalueTolerances(iTol) ...
            * max(incrementSpectralNorm(k), referenceIncrementScale);

        positiveRank(k, iTol) = ...
            nnz(lambda > threshold);
    end

    constructionThreshold = ...
        positivePartTolerance ...
        * max(incrementSpectralNorm(k), referenceIncrementScale);

    lambdaPositive = max(lambda, 0);
    lambdaNegative = max(-lambda, 0);

    % Remove positive and negative values below the construction threshold.
    lambdaPositive(lambdaPositive <= constructionThreshold) = 0;
    lambdaNegative(lambdaNegative <= constructionThreshold) = 0;

    retainedPositiveRank(k) = nnz(lambdaPositive > 0);

    Dplus = V * diag(lambdaPositive) * V.';
    Dminus = V * diag(lambdaNegative) * V.';

    Dplus = 0.5 * (Dplus + Dplus.');
    Dminus = 0.5 * (Dminus + Dminus.');

    positivePart{k} = Dplus;
    negativePart{k} = Dminus;

    positiveTrace(k) = trace(Dplus);
    negativeTrace(k) = trace(Dminus);

    positiveSpectralNorm(k) = max(lambdaPositive);
    negativeSpectralNorm(k) = max(lambdaNegative);

    positiveFrobeniusNorm(k) = norm(Dplus, 'fro');
    negativeFrobeniusNorm(k) = norm(Dminus, 'fro');
end

% -------------------------------------------------------------------------
% Construct cumulative positive-growth matrix.
% -------------------------------------------------------------------------

Kpositive = zeros(r, r);
aggregateIntervalWeight = zeros(numIntervals, 1);

for k = 1:numIntervals
    Dplus = positivePart{k};

    switch aggregateWeightMethod
        case 'magnitude'
            % If D is per unit x, multiply by dx so the aggregate measures
            % total positive variation over each interval.
            if strcmp(incrementScaling, 'perunitx')
                intervalWeight = dx(k);
            else
                intervalWeight = 1;
            end

        case 'spectralnormalized'
            if positiveSpectralNorm(k) > 0
                intervalWeight = 1 / positiveSpectralNorm(k);
            else
                intervalWeight = 0;
            end

        case 'tracenormalized'
            if positiveTrace(k) > 0
                intervalWeight = 1 / positiveTrace(k);
            else
                intervalWeight = 0;
            end

        otherwise
            error('analyzePositiveIncrementSubspace:UnknownWeightMethod', ...
                'Unknown aggregateWeightMethod "%s".', ...
                options.aggregateWeightMethod);
    end

    aggregateIntervalWeight(k) = intervalWeight;

    Kpositive = ...
        Kpositive + intervalWeight * Dplus;
end

Kpositive = 0.5 * (Kpositive + Kpositive.');

[Vgrowth, Dgrowth] = eig(Kpositive);

growthEigenvalues = real(diag(Dgrowth));

[growthEigenvalues, order] = ...
    sort(growthEigenvalues, 'descend');

Vgrowth = real(Vgrowth(:, order));

positiveGrowthEigenvalues = max(growthEigenvalues, 0);
totalGrowthTrace = sum(positiveGrowthEigenvalues);

if totalGrowthTrace > 0
    cumulativeGrowthFraction = ...
        cumsum(positiveGrowthEigenvalues) / totalGrowthTrace;
else
    cumulativeGrowthFraction = ...
        zeros(size(positiveGrowthEigenvalues));
end

% -------------------------------------------------------------------------
% Cumulative numerical dimension of all positive increment ranges.
% -------------------------------------------------------------------------

cumulativePositiveMatrix = zeros(r, r);

cumulativePositiveEigenvalues = zeros(r, numIntervals);
cumulativePositiveRank = zeros( ...
    numIntervals, numel(positiveEigenvalueTolerances));

for k = 1:numIntervals
    Dplus = positivePart{k};

    scale = max(positiveSpectralNorm(k), eps);

    % Normalize each nonzero interval contribution so that intervals with
    % weak growth still contribute their range to the span diagnostic.
    if scale > 0
        cumulativePositiveMatrix = ...
            cumulativePositiveMatrix + Dplus / scale;
    end

    cumulativePositiveMatrix = ...
        0.5 * ( ...
            cumulativePositiveMatrix ...
            + cumulativePositiveMatrix.');

    lambdaCumulative = sort( ...
        real(eig(cumulativePositiveMatrix)), 'descend');

    cumulativePositiveEigenvalues(:, k) = ...
        lambdaCumulative;

    for iTol = 1:numel(positiveEigenvalueTolerances)
        cumulativePositiveRank(k, iTol) = ...
            relativeRank( ...
                lambdaCumulative, ...
                positiveEigenvalueTolerances(iTol));
    end
end

% -------------------------------------------------------------------------
% Candidate fixed-basis diagnostics.
% -------------------------------------------------------------------------

numRanks = numel(candidateRanks);

aggregateCapturedFraction = zeros(numRanks, 1);

maximumPositiveOperatorResidual = zeros(numRanks, 1);
medianPositiveOperatorResidual = zeros(numRanks, 1);

weightedPositiveTraceCapture = zeros(numRanks, 1);
minimumIntervalPositiveTraceCapture = zeros(numRanks, 1);

maximumComplementGrowth = zeros(numRanks, 1);
maximumRelativeComplementGrowth = zeros(numRanks, 1);
numComplementGrowthViolations = zeros(numRanks, 1);

maximumCrossCoupling = zeros(numRanks, 1);
maximumRelativeCrossCoupling = zeros(numRanks, 1);

complementGrowthByInterval = zeros(numIntervals, numRanks);
relativeComplementGrowthByInterval = zeros(numIntervals, numRanks);

crossCouplingByInterval = zeros(numIntervals, numRanks);
relativeCrossCouplingByInterval = zeros(numIntervals, numRanks);

positiveResidualByInterval = zeros(numIntervals, numRanks);
positiveTraceCaptureByInterval = ones(numIntervals, numRanks);

totalWeightedPositiveTrace = sum( ...
    aggregateIntervalWeight .* positiveTrace);

for iRank = 1:numRanks
    p = candidateRanks(iRank);

    P = Vgrowth(:, 1:p);
    projector = P * P.';
    complementProjector = eye(r) - projector;

    aggregateCapturedFraction(iRank) = ...
        sum(positiveGrowthEigenvalues(1:p)) ...
        / max(totalGrowthTrace, eps);

    weightedCapturedTrace = 0;

    for k = 1:numIntervals
        Dk = increment{k};
        Dplus = positivePart{k};

        localIncrementScale = max( ...
            incrementSpectralNorm(k), eps);

        localPositiveScale = max( ...
            positiveSpectralNorm(k), eps);

        DplusProjected = ...
            projector * Dplus * projector;

        positiveResidualByInterval(k, iRank) = ...
            norm(Dplus - DplusProjected, 2) ...
            / localPositiveScale;

        if positiveTrace(k) > 0
            capturedTrace = trace(P.' * Dplus * P);

            positiveTraceCaptureByInterval(k, iRank) = ...
                capturedTrace / positiveTrace(k);

            weightedCapturedTrace = ...
                weightedCapturedTrace ...
                + aggregateIntervalWeight(k) * capturedTrace;
        else
            positiveTraceCaptureByInterval(k, iRank) = 1;
        end

        complementIncrement = ...
            complementProjector * Dk * complementProjector;

        complementIncrement = ...
            0.5 * ( ...
                complementIncrement ...
                + complementIncrement.');

        lambdaComplementMaximum = max(real(eig( ...
            complementIncrement)));

        complementGrowthByInterval(k, iRank) = ...
            lambdaComplementMaximum;

        relativeComplementGrowthByInterval(k, iRank) = ...
            lambdaComplementMaximum / localIncrementScale;

        crossBlock = ...
            P.' * Dk * complementProjector;

        crossCouplingByInterval(k, iRank) = ...
            norm(crossBlock, 2);

        relativeCrossCouplingByInterval(k, iRank) = ...
            norm(crossBlock, 2) / localIncrementScale;
    end

    maximumPositiveOperatorResidual(iRank) = ...
        max(positiveResidualByInterval(:, iRank));

    medianPositiveOperatorResidual(iRank) = ...
        median(positiveResidualByInterval(:, iRank));

    weightedPositiveTraceCapture(iRank) = ...
        weightedCapturedTrace ...
        / max(totalWeightedPositiveTrace, eps);

    intervalsWithPositiveTrace = positiveTrace > 0;

    if any(intervalsWithPositiveTrace)
        minimumIntervalPositiveTraceCapture(iRank) = ...
            min(positiveTraceCaptureByInterval( ...
                intervalsWithPositiveTrace, iRank));
    else
        minimumIntervalPositiveTraceCapture(iRank) = 1;
    end

    maximumComplementGrowth(iRank) = ...
        max(complementGrowthByInterval(:, iRank));

    maximumRelativeComplementGrowth(iRank) = ...
        max(relativeComplementGrowthByInterval(:, iRank));

    numComplementGrowthViolations(iRank) = nnz( ...
        complementGrowthByInterval(:, iRank) ...
        > complementTolerance * referenceIncrementScale);

    maximumCrossCoupling(iRank) = ...
        max(crossCouplingByInterval(:, iRank));

    maximumRelativeCrossCoupling(iRank) = ...
        max(relativeCrossCouplingByInterval(:, iRank));
end

numericallyMonotoneComplement = ...
    numComplementGrowthViolations == 0;

% -------------------------------------------------------------------------
% Package report.
% -------------------------------------------------------------------------

report = struct();

report.eventType = 'discreteAdjacentIncrementDiagnostic';
report.continuousPathCertified = false;
report.exactCertificateAvailable = false;

report.x = x;
report.intervalCoordinate = ...
    0.5 * (x(1:end - 1) + x(2:end));
report.dx = dx;

report.fullStochasticDimension = r;
report.numStations = Nx;
report.numIntervals = numIntervals;

report.incrementScaling = incrementScaling;
report.increment = increment;
report.incrementEigenvalues = incrementEigenvalues;
report.incrementSpectralNorm = incrementSpectralNorm;
report.incrementFrobeniusNorm = incrementFrobeniusNorm;

report.positiveEigenvalueTolerances = ...
    positiveEigenvalueTolerances;

report.positivePartTolerance = ...
    positivePartTolerance;

report.positivePart = positivePart;
report.negativePart = negativePart;

report.positiveRank = positiveRank;
report.retainedPositiveRank = retainedPositiveRank;

report.positiveTrace = positiveTrace;
report.negativeTrace = negativeTrace;

report.positiveSpectralNorm = positiveSpectralNorm;
report.negativeSpectralNorm = negativeSpectralNorm;

report.positiveFrobeniusNorm = positiveFrobeniusNorm;
report.negativeFrobeniusNorm = negativeFrobeniusNorm;

report.aggregateWeightMethod = aggregateWeightMethod;
report.aggregateIntervalWeight = aggregateIntervalWeight;

report.aggregatePositiveMatrix = Kpositive;
report.growthBasis = Vgrowth;
report.growthEigenvalues = growthEigenvalues;
report.cumulativeGrowthFraction = cumulativeGrowthFraction;

report.cumulativePositiveEigenvalues = ...
    cumulativePositiveEigenvalues;
report.cumulativePositiveRank = ...
    cumulativePositiveRank;

report.candidateRanks = candidateRanks;
report.aggregateCapturedFraction = ...
    aggregateCapturedFraction;

report.maximumPositiveOperatorResidual = ...
    maximumPositiveOperatorResidual;
report.medianPositiveOperatorResidual = ...
    medianPositiveOperatorResidual;

report.weightedPositiveTraceCapture = ...
    weightedPositiveTraceCapture;
report.minimumIntervalPositiveTraceCapture = ...
    minimumIntervalPositiveTraceCapture;

report.complementTolerance = complementTolerance;
report.referenceIncrementScale = ...
    referenceIncrementScale;

report.complementGrowthByInterval = ...
    complementGrowthByInterval;
report.relativeComplementGrowthByInterval = ...
    relativeComplementGrowthByInterval;

report.maximumComplementGrowth = ...
    maximumComplementGrowth;
report.maximumRelativeComplementGrowth = ...
    maximumRelativeComplementGrowth;

report.numComplementGrowthViolations = ...
    numComplementGrowthViolations;
report.numericallyMonotoneComplement = ...
    numericallyMonotoneComplement;

report.crossCouplingByInterval = ...
    crossCouplingByInterval;
report.relativeCrossCouplingByInterval = ...
    relativeCrossCouplingByInterval;

report.maximumCrossCoupling = ...
    maximumCrossCoupling;
report.maximumRelativeCrossCoupling = ...
    maximumRelativeCrossCoupling;

report.positiveResidualByInterval = ...
    positiveResidualByInterval;
report.positiveTraceCaptureByInterval = ...
    positiveTraceCaptureByInterval;

report.interpretation = ...
    ['A fixed basis is constructed from positive adjacent increments ', ...
     'of the real stochastic energy matrices. Complement monotonicity ', ...
     'is tested numerically at stored streamwise intervals only.'];

if verbose
    fprintf('\nPositive-increment growth-subspace analysis:\n');
    fprintf('  full stochastic dimension = %d\n', r);
    fprintf('  number of intervals       = %d\n', numIntervals);
    fprintf('  increment scaling         = %s\n', incrementScaling);
    fprintf('  aggregate weighting       = %s\n', ...
        aggregateWeightMethod);
    fprintf('  positive-part tolerance   = %.3e\n', ...
        positivePartTolerance);

    fprintf('\n  Positive eigenvalue counts:\n');

    for iTol = 1:numel(positiveEigenvalueTolerances)
        ranks = positiveRank(:, iTol);

        fprintf('    tolerance %.1e: min=%d, median=%.1f, max=%d\n', ...
            positiveEigenvalueTolerances(iTol), ...
            min(ranks), median(ranks), max(ranks));
    end

    fprintf('\n  Final cumulative positive-range ranks:\n');

    for iTol = 1:numel(positiveEigenvalueTolerances)
        fprintf('    tolerance %.1e: %d\n', ...
            positiveEigenvalueTolerances(iTol), ...
            cumulativePositiveRank(end, iTol));
    end

    fprintf('\n  Candidate growth bases:\n');
    fprintf([ ...
        '    %6s %12s %12s %12s %12s %12s %10s\n'], ...
        'rank', ...
        'agg frac', ...
        'trace cap', ...
        'max +resid', ...
        'max Qgrowth', ...
        'max coupling', ...
        'violations');

    for iRank = 1:numRanks
        fprintf([ ...
            '    %6d %12.5e %12.5e %12.5e ', ...
            '%12.5e %12.5e %10d\n'], ...
            candidateRanks(iRank), ...
            aggregateCapturedFraction(iRank), ...
            weightedPositiveTraceCapture(iRank), ...
            maximumPositiveOperatorResidual(iRank), ...
            maximumRelativeComplementGrowth(iRank), ...
            maximumRelativeCrossCoupling(iRank), ...
            numComplementGrowthViolations(iRank));
    end

    monotoneIndices = find(numericallyMonotoneComplement);

    if isempty(monotoneIndices)
        fprintf('\n  No tested reduced rank produced a numerically ');
        fprintf('monotone complement.\n');
    else
        fprintf('\n  Smallest tested rank with numerically monotone ');
        fprintf('complement: %d\n', ...
            candidateRanks(monotoneIndices(1)));
    end

    fprintf('\n');
end

end

% =========================================================================
function rankValue = relativeRank(eigenvalues, tolerance)

eigenvalues = real(eigenvalues(:));

scale = max(abs(eigenvalues));

if scale == 0
    rankValue = 0;
else
    rankValue = nnz( ...
        eigenvalues > tolerance * scale);
end

end

% =========================================================================
function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end