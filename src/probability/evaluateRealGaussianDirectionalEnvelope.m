function result = evaluateRealGaussianDirectionalEnvelope( ...
    gain, weights, stochasticDimension, threshold)
%EVALUATEREALGAUSSIANDIRECTIONALENVELOPE Evaluate local and first passage.
%
% Inputs:
%   gain                 K-by-Nx directional quadratic forms
%   weights              K-by-1 angular weights
%   stochasticDimension  chi-square degrees of freedom
%   threshold            positive scalar constant threshold
%
% Output:
%   localProbability
%   runningMaximum
%   transitionCdf
%   survival
%
% Negative instantaneous gains are permitted because an approximate matrix
% family may be indefinite. The first-passage running maximum is initialized
% from the first station and remains usable if it is positive.
%
% MATLAB version: R2020b

[K, Nx] = size(gain);
weights = weights(:);

if numel(weights) ~= K
    error('evaluateRealGaussianDirectionalEnvelope:WeightMismatch', ...
        'The angular weights are incompatible with gain.');
end

validateattributes(stochasticDimension, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'stochasticDimension');

validateattributes(threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'threshold');

runningMaximum = zeros(K, Nx);
runningMaximum(:, 1) = gain(:, 1);

for n = 2:Nx
    runningMaximum(:, n) = max( ...
        runningMaximum(:, n - 1), ...
        gain(:, n));
end

localProbability = zeros(Nx, 1);
transitionCdf = zeros(Nx, 1);

for n = 1:Nx
    localGain = gain(:, n);
    historyGain = runningMaximum(:, n);

    localPositive = localGain > 0;
    historyPositive = historyGain > 0;

    localTail = zeros(K, 1);
    historyTail = zeros(K, 1);

    localTail(localPositive) = gammainc( ...
        threshold ./ (2 * localGain(localPositive)), ...
        stochasticDimension / 2, ...
        'upper');

    historyTail(historyPositive) = gammainc( ...
        threshold ./ (2 * historyGain(historyPositive)), ...
        stochasticDimension / 2, ...
        'upper');

    localProbability(n) = weights.' * localTail;
    transitionCdf(n) = weights.' * historyTail;
end

result = struct();

result.gain = gain;
result.runningMaximum = runningMaximum;

result.localProbability = localProbability;
result.transitionCdf = transitionCdf;
result.survival = 1 - transitionCdf;

result.pCensored = result.survival(end);

end