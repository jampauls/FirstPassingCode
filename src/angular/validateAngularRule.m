function diagnostics = validateAngularRule(U, weights)
%VALIDATEANGULARRULE Validate sphere directions and angular weights.

if ~isnumeric(U) || ndims(U) ~= 2 || isempty(U)
    error('validateAngularRule:InvalidDirections', ...
        'U must be a nonempty numeric matrix.');
end

weights = weights(:);
K = size(U, 2);

if numel(weights) ~= K
    error('validateAngularRule:WeightLengthMismatch', ...
        'The number of weights must equal the number of directions.');
end

if any(~isfinite(U(:))) || any(~isfinite(weights))
    error('validateAngularRule:NonfiniteValues', ...
        'Angular directions and weights must be finite.');
end

norms = sqrt(sum(abs(U).^2, 1));
normDefect = max(abs(norms - 1));

if normDefect > 1e-10
    error('validateAngularRule:NonunitDirections', ...
        'Maximum direction norm defect is %.3e.', normDefect);
end

if any(weights < 0)
    error('validateAngularRule:NegativeWeights', ...
        'Angular weights must be nonnegative.');
end

weightDefect = abs(sum(weights) - 1);

if weightDefect > 1e-12
    error('validateAngularRule:WeightNormalization', ...
        'Angular weights sum to %.16g instead of one.', sum(weights));
end

% A useful isotropy diagnostic:
% E[u*u'] = I/r for a uniform sphere.
r = size(U, 1);
secondMoment = U * diag(weights) * U';
secondMomentTarget = eye(r) / r;

isotropyDefect = norm( ...
    secondMoment - secondMomentTarget, 'fro');

diagnostics = struct();
diagnostics.maxNormDefect = normDefect;
diagnostics.weightDefect = weightDefect;
diagnostics.secondMoment = secondMoment;
diagnostics.isotropyDefect = isotropyDefect;

end
