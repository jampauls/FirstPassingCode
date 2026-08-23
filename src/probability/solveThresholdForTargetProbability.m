function solution = solveThresholdForTargetProbability( ...
    envelope, targetProbability, amplitude, options)
%SOLVETHRESHOLDFORTARGETPROBABILITY Solve a constant threshold.
%
% Finds threshold e such that:
%
%   F_Xtr(xTarget; e, amplitude) = targetProbability.
%
% By default, xTarget is the final streamwise station.
%
% Since transition probability decreases monotonically with threshold,
% logarithmic bisection is used.
%
% Options:
%   targetIndex
%   initialThreshold
%   lowerThreshold
%   upperThreshold
%   probabilityTolerance
%   relativeTolerance
%   maxIterations
%   verbose
%
% MATLAB version: R2020b

if nargin < 4
    options = struct();
end

validateattributes(targetProbability, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1}, ...
    mfilename, 'targetProbability');

validateattributes(amplitude, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'amplitude');

targetIndex = getOption( ...
    options, 'targetIndex', envelope.numStations);

validateattributes(targetIndex, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, ...
     '<=', envelope.numStations}, ...
    mfilename, 'options.targetIndex');

probabilityTolerance = getOption( ...
    options, 'probabilityTolerance', 1e-6);

relativeTolerance = getOption( ...
    options, 'relativeTolerance', 1e-6);

maxIterations = getOption( ...
    options, 'maxIterations', 100);

verbose = getOption( ...
    options, 'verbose', true);

initialThreshold = getOption( ...
    options, 'initialThreshold', []);

if isempty(initialThreshold)
    % Construct a characteristic energy scale from the terminal running
    % maxima over all directions and replicates.
    characteristicValues = zeros( ...
        envelope.numReplicates, 1);

    for j = 1:envelope.numReplicates
        characteristicValues(j) = median( ...
            envelope.mReplicate{j}(:, targetIndex));
    end

    initialThreshold = ...
        amplitude^2 * median(characteristicValues);

    if initialThreshold <= 0 || ~isfinite(initialThreshold)
        initialThreshold = 1;
    end
end

lowerThreshold = getOption( ...
    options, 'lowerThreshold', ...
    initialThreshold / 10);

upperThreshold = getOption( ...
    options, 'upperThreshold', ...
    initialThreshold * 10);

% -------------------------------------------------------------------------
% Expand the bracket if necessary.
% -------------------------------------------------------------------------

[probabilityLower, seLower] = ...
    evaluateTargetProbability( ...
        envelope, lowerThreshold, ...
        amplitude, targetIndex);

[probabilityUpper, seUpper] = ...
    evaluateTargetProbability( ...
        envelope, upperThreshold, ...
        amplitude, targetIndex);

numExpansion = 0;
maximumExpansion = 50;

while probabilityLower < targetProbability && ...
        numExpansion < maximumExpansion

    lowerThreshold = lowerThreshold / 10;

    [probabilityLower, seLower] = ...
        evaluateTargetProbability( ...
            envelope, lowerThreshold, ...
            amplitude, targetIndex);

    numExpansion = numExpansion + 1;
end

while probabilityUpper > targetProbability && ...
        numExpansion < maximumExpansion

    upperThreshold = upperThreshold * 10;

    [probabilityUpper, seUpper] = ...
        evaluateTargetProbability( ...
            envelope, upperThreshold, ...
            amplitude, targetIndex);

    numExpansion = numExpansion + 1;
end

if probabilityLower < targetProbability || ...
        probabilityUpper > targetProbability
    error('solveThresholdForTargetProbability:BracketFailure', ...
        ['Unable to bracket the requested target probability. ', ...
         'F(lower)=%.6e and F(upper)=%.6e.'], ...
        probabilityLower, probabilityUpper);
end

% -------------------------------------------------------------------------
% Logarithmic bisection
% -------------------------------------------------------------------------

historyThreshold = zeros(maxIterations, 1);
historyProbability = zeros(maxIterations, 1);
historyStandardError = zeros(maxIterations, 1);

converged = false;

for iteration = 1:maxIterations
    thresholdMid = sqrt( ...
        lowerThreshold * upperThreshold);

    [probabilityMid, seMid] = ...
        evaluateTargetProbability( ...
            envelope, thresholdMid, ...
            amplitude, targetIndex);

    historyThreshold(iteration) = thresholdMid;
    historyProbability(iteration) = probabilityMid;
    historyStandardError(iteration) = seMid;

    probabilityError = ...
        probabilityMid - targetProbability;

    relativeBracketWidth = ...
        upperThreshold / lowerThreshold - 1;

    if abs(probabilityError) <= probabilityTolerance || ...
            relativeBracketWidth <= relativeTolerance
        converged = true;
        break;
    end

    % Increasing the threshold decreases transition probability.
    if probabilityMid > targetProbability
        lowerThreshold = thresholdMid;
        probabilityLower = probabilityMid;
        seLower = seMid;
    else
        upperThreshold = thresholdMid;
        probabilityUpper = probabilityMid;
        seUpper = seMid;
    end
end

historyThreshold = historyThreshold(1:iteration);
historyProbability = historyProbability(1:iteration);
historyStandardError = historyStandardError(1:iteration);

finalResult = evaluateEnergyEnvelope( ...
    envelope, thresholdMid, amplitude);

solution = struct();

solution.threshold = thresholdMid;
solution.amplitude = amplitude;
solution.targetProbability = targetProbability;
solution.targetIndex = targetIndex;
solution.targetCoordinate = envelope.x(targetIndex);

solution.achievedProbability = ...
    finalResult.F(targetIndex);

solution.standardError = ...
    finalResult.standardError(targetIndex);

solution.probabilityError = ...
    solution.achievedProbability - targetProbability;

solution.converged = converged;
solution.iterations = iteration;
solution.lowerThreshold = lowerThreshold;
solution.upperThreshold = upperThreshold;

solution.historyThreshold = historyThreshold;
solution.historyProbability = historyProbability;
solution.historyStandardError = historyStandardError;

solution.probabilityResult = finalResult;

if verbose
    fprintf('\nConstant-threshold solution:\n');
    fprintf('  target coordinate      = %.8e\n', ...
        solution.targetCoordinate);
    fprintf('  target probability     = %.8e\n', ...
        targetProbability);
    fprintf('  amplitude multiplier   = %.8e\n', ...
        amplitude);
    fprintf('  threshold              = %.8e\n', ...
        solution.threshold);
    fprintf('  achieved probability   = %.8e\n', ...
        solution.achievedProbability);
    fprintf('  RQMC standard error    = %.8e\n', ...
        solution.standardError);
    fprintf('  probability residual   = %.8e\n', ...
        solution.probabilityError);
    fprintf('  iterations             = %d\n', ...
        solution.iterations);
    fprintf('  converged              = %d\n\n', ...
        solution.converged);
end

end

% =========================================================================
function [probability, standardError] = ...
    evaluateTargetProbability( ...
        envelope, threshold, amplitude, targetIndex)

r = envelope.r;
J = envelope.numReplicates;
weights = envelope.weights(:);

replicateProbability = zeros(J, 1);
amplitudeSquared = amplitude^2;

for j = 1:J
    m = envelope.mReplicate{j}(:, targetIndex);

    positive = m > 0;
    directionalProbability = zeros(size(m));

    directionalProbability(positive) = ...
        gammainc( ...
            threshold ...
            ./ (2 * amplitudeSquared * m(positive)), ...
            r / 2, ...
            'upper');

    replicateProbability(j) = ...
        weights.' * directionalProbability;
end

probability = mean(replicateProbability);
standardError = std(replicateProbability) / sqrt(J);

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