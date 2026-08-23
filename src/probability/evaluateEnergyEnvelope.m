function result = evaluateEnergyEnvelope( ...
    envelope, threshold, amplitude)
%EVALUATEENERGYENVELOPE Evaluate one threshold/amplitude combination.
%
%   result = evaluateEnergyEnvelope(envelope, threshold, amplitude)
%
% Transition condition:
%
%   amplitude^2 * R^2 * M_u(x) >= threshold.
%
% Therefore:
%
%   P(Xtr <= x | u)
%       = chi2sf(threshold / (amplitude^2*M_u(x)); r).
%
% The local exceedance probability is evaluated using the endpoint gains
% stored in the same RQMC replicates.
%
% MATLAB version: R2020b

validateattributes(threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'threshold');

validateattributes(amplitude, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'amplitude');

r = envelope.r;
J = envelope.numReplicates;
Nx = envelope.numStations;
weights = envelope.weights(:);

amplitudeSquared = amplitude^2;

FReplicate = zeros(Nx, J);
pLocalReplicate = zeros(Nx, J);

for j = 1:J
    m = envelope.mReplicate{j};
    a = envelope.aReplicate{j};

    FReplicate(:, j) = evaluateDirectionalTail( ...
        m, weights, r, threshold, amplitudeSquared);

    pLocalReplicate(:, j) = evaluateDirectionalTail( ...
        a, weights, r, threshold, amplitudeSquared);
end

F = mean(FReplicate, 2);
S = 1 - F;
pLocal = mean(pLocalReplicate, 2);

memoryCorrectionReplicate = ...
    FReplicate - pLocalReplicate;

memoryCorrection = ...
    mean(memoryCorrectionReplicate, 2);

standardError = ...
    std(FReplicate, 0, 2) / sqrt(J);

pLocalStandardError = ...
    std(pLocalReplicate, 0, 2) / sqrt(J);

memoryCorrectionStandardError = ...
    std(memoryCorrectionReplicate, 0, 2) / sqrt(J);

pInterval = zeros(Nx, 1);
pInterval(1) = F(1);

if Nx > 1
    pInterval(2:end) = diff(F);
end

smallNegative = ...
    pInterval < 0 & pInterval > -1e-13;

pInterval(smallNegative) = 0;

if any(pInterval < -1e-12)
    error('evaluateEnergyEnvelope:NonmonotoneCDF', ...
        'The replicate-mean transition CDF is materially nonmonotone.');
end

localViolation = max(pLocal - F);

if localViolation > 1e-12
    error('evaluateEnergyEnvelope:LocalProbabilityViolation', ...
        ['The local exceedance probability exceeds the first-transition ', ...
         'CDF by %.3e.'], localViolation);
end

result = struct();

result.threshold = threshold;
result.amplitude = amplitude;
result.effectiveThreshold = ...
    threshold / amplitudeSquared;

result.F = F;
result.S = S;
result.pLocal = pLocal;
result.memoryCorrection = memoryCorrection;

result.standardError = standardError;
result.pLocalStandardError = ...
    pLocalStandardError;
result.memoryCorrectionStandardError = ...
    memoryCorrectionStandardError;

result.FReplicate = FReplicate;
result.pLocalReplicate = pLocalReplicate;
result.memoryCorrectionReplicate = ...
    memoryCorrectionReplicate;

result.pInterval = pInterval;
result.pCensored = S(end);
result.terminalStandardError = ...
    standardError(end);

end

% =========================================================================
function probability = evaluateDirectionalTail( ...
    gain, weights, r, threshold, amplitudeSquared)
%EVALUATEDIRECTIONALTAIL Evaluate angular chi-square tail averages.
%
% gain is K-by-Nx.

[K, Nx] = size(gain);

if numel(weights) ~= K
    error('evaluateDirectionalTail:WeightLengthMismatch', ...
        'The weight vector is incompatible with the gain array.');
end

if any(~isfinite(gain(:)))
    error('evaluateDirectionalTail:NonfiniteGain', ...
        'The directional gains contain nonfinite values.');
end

gainScale = max(max(abs(gain(:))), 1);

smallNegative = ...
    gain < 0 & gain >= -1e-12 * gainScale;

gain(smallNegative) = 0;

if any(gain(:) < 0)
    error('evaluateDirectionalTail:NegativeGain', ...
        'The directional gain contains material negative values.');
end

probability = zeros(Nx, 1);

for n = 1:Nx
    gn = gain(:, n);
    positive = gn > 0;

    directionalProbability = zeros(K, 1);

    argument = ...
        threshold ...
        ./ (2 * amplitudeSquared * gn(positive));

    directionalProbability(positive) = ...
        gammainc(argument, r / 2, 'upper');

    probability(n) = ...
        weights.' * directionalProbability;
end

probability = min(max(probability, 0), 1);

end