function pLocal = computeLocalExceedanceFromGains(a, weights, r)
%COMPUTELOCALEXCEEDANCEFROMGAINS Compute local exceedance probabilities.
%
%   pLocal = computeLocalExceedanceFromGains(a, weights, r)
%
% Inputs:
%   a        K-by-Nx directional endpoint gains:
%
%                a(k,n) = u_k' * G(x_n) * u_k
%
%   weights  K-by-1 normalized angular weights
%   r        real Gaussian stochastic dimension
%
% Output:
%   pLocal   Nx-by-1 local threshold-exceedance probability
%
% For each direction:
%
%   P(R^2 * a(k,n) >= 1)
%       = chi2sf(1/a(k,n); r).
%
% MATLAB version: R2020b

if ~isnumeric(a) || ndims(a) ~= 2 || isempty(a)
    error('computeLocalExceedanceFromGains:InvalidGainArray', ...
        'a must be a nonempty K-by-Nx numeric matrix.');
end

[K, Nx] = size(a);
weights = weights(:);

if numel(weights) ~= K
    error('computeLocalExceedanceFromGains:WeightLengthMismatch', ...
        'numel(weights) must equal size(a,1).');
end

if any(weights < 0) || ...
        abs(sum(weights) - 1) > 1e-12
    error('computeLocalExceedanceFromGains:InvalidWeights', ...
        'Weights must be nonnegative and sum to one.');
end

validateattributes(r, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'r');

if any(~isfinite(a(:)))
    error('computeLocalExceedanceFromGains:NonfiniteGain', ...
        'The directional gains contain nonfinite values.');
end

% Remove only roundoff-scale negative values.
gainScale = max(max(abs(a(:))), 1);
smallNegative = ...
    a < 0 & a >= -1e-12 * gainScale;

a(smallNegative) = 0;

if any(a(:) < 0)
    error('computeLocalExceedanceFromGains:NegativeGain', ...
        'The directional gains contain material negative values.');
end

pLocal = zeros(Nx, 1);

for n = 1:Nx
    gain = a(:, n);
    positive = gain > 0;

    directionalProbability = zeros(K, 1);

    directionalProbability(positive) = ...
        gammainc( ...
            1 ./ (2 * gain(positive)), ...
            r / 2, ...
            'upper');

    pLocal(n) = ...
        weights.' * directionalProbability;
end

pLocal = min(max(pLocal, 0), 1);

end
