function diagnostics = computeDirectionalDiagnostics(m, U, weights, r, cfg)
%COMPUTEDIRECTIONALDIAGNOSTICS Summarize angular transition contributions.
%
% Inputs:
%   m        K-by-Nx running maxima
%   U        r-by-K direction matrix
%   weights  K-by-1 angular weights
%
% Output includes:
%   terminalConditionalProbability
%   terminalWeightedContribution
%   transitionConditionedWeights
%   transitionDirectionSecondMoment
%   dominantContributionIndices

if nargin < 5
    cfg = struct(); %#ok<NASGU>
end

[K, Nx] = size(m);
weights = weights(:);

if size(U, 2) ~= K || size(U, 1) ~= r
    error('computeDirectionalDiagnostics:DimensionMismatch', ...
        'Dimensions of m and U are inconsistent.');
end

mTerminal = m(:, Nx);

positive = mTerminal > 0;
psi = zeros(K, 1);

psi(positive) = ...
    gammainc( ...
        1 ./ (2 * mTerminal(positive)), ...
        r / 2, ...
        'upper');

contribution = weights .* psi;
terminalCdf = sum(contribution);

if terminalCdf > 0
    conditionalWeights = contribution / terminalCdf;
else
    conditionalWeights = zeros(K, 1);
end

secondMoment = zeros(r, r);

for k = 1:K
    secondMoment = secondMoment ...
        + conditionalWeights(k) * (U(:, k) * U(:, k).');
end

[sortedContribution, order] = sort(contribution, 'descend');

numReported = min(10, K);

diagnostics = struct();
diagnostics.terminalRunningMaximum = mTerminal;
diagnostics.terminalConditionalProbability = psi;
diagnostics.terminalWeightedContribution = contribution;
diagnostics.transitionConditionedWeights = conditionalWeights;
diagnostics.transitionDirectionSecondMoment = secondMoment;
diagnostics.terminalCdf = terminalCdf;

diagnostics.dominantContributionIndices = order(1:numReported);
diagnostics.dominantContributions = sortedContribution(1:numReported);
diagnostics.dominantDirections = U(:, order(1:numReported));

end
