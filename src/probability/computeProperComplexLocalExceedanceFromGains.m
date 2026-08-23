function pLocal = ...
    computeProperComplexLocalExceedanceFromGains( ...
        gain, weights, complexDimension)
%COMPUTEPROPERCOMPLEXLOCALEXCEEDANCEFROMGAINS Local proper-complex tails.

[K, Nx] = size(gain);
weights = weights(:);

if numel(weights) ~= K
    error('computeProperComplexLocalExceedanceFromGains:WeightMismatch', ...
        'The weight vector is incompatible with the gain array.');
end

gainScale = max(max(abs(gain(:))), 1);

smallNegative = ...
    gain < 0 & gain >= -1e-12 * gainScale;

gain(smallNegative) = 0;

if any(gain(:) < 0)
    error('computeProperComplexLocalExceedanceFromGains:NegativeGain', ...
        'A material negative directional gain was detected.');
end

pLocal = zeros(Nx, 1);

for n = 1:Nx
    gn = gain(:, n);
    positive = gn > 0;

    directionalProbability = zeros(K, 1);

    directionalProbability(positive) = ...
        gammainc( ...
            1 ./ gn(positive), ...
            complexDimension, ...
            'upper');

    pLocal(n) = ...
        weights.' * directionalProbability;
end

pLocal = min(max(pLocal, 0), 1);

end
