function family = evaluateThresholdAmplitudeFamily( ...
    envelope, thresholdValues, amplitudeValues)
%EVALUATETHRESHOLDAMPLITUDEFAMILY Evaluate a parameter family.
%
% Outputs are stored as:
%
%   F(:,iThreshold,iAmplitude)
%
% MATLAB version: R2020b

thresholdValues = thresholdValues(:);
amplitudeValues = amplitudeValues(:);

if isempty(thresholdValues) || ...
        any(~isfinite(thresholdValues)) || ...
        any(thresholdValues <= 0)
    error('evaluateThresholdAmplitudeFamily:InvalidThresholds', ...
        'Threshold values must be finite and positive.');
end

if isempty(amplitudeValues) || ...
        any(~isfinite(amplitudeValues)) || ...
        any(amplitudeValues <= 0)
    error('evaluateThresholdAmplitudeFamily:InvalidAmplitudes', ...
        'Amplitude values must be finite and positive.');
end

Nx = envelope.numStations;
nThreshold = numel(thresholdValues);
nAmplitude = numel(amplitudeValues);

F = zeros(Nx, nThreshold, nAmplitude);
S = zeros(Nx, nThreshold, nAmplitude);
pLocal = zeros(Nx, nThreshold, nAmplitude);
memoryCorrection = zeros(Nx, nThreshold, nAmplitude);

standardError = zeros(Nx, nThreshold, nAmplitude);
memoryCorrectionStandardError = ...
    zeros(Nx, nThreshold, nAmplitude);

pCensored = zeros(nThreshold, nAmplitude);
terminalProbability = zeros(nThreshold, nAmplitude);
terminalStandardError = zeros(nThreshold, nAmplitude);

fprintf('\nEvaluating threshold/amplitude family:\n');
fprintf('  thresholds = %d\n', nThreshold);
fprintf('  amplitudes = %d\n\n', nAmplitude);

for iAmplitude = 1:nAmplitude
    for iThreshold = 1:nThreshold
        result = evaluateEnergyEnvelope( ...
            envelope, ...
            thresholdValues(iThreshold), ...
            amplitudeValues(iAmplitude));

        F(:, iThreshold, iAmplitude) = result.F;
        S(:, iThreshold, iAmplitude) = result.S;
        pLocal(:, iThreshold, iAmplitude) = result.pLocal;

        memoryCorrection(:, iThreshold, iAmplitude) = ...
            result.memoryCorrection;

        standardError(:, iThreshold, iAmplitude) = ...
            result.standardError;

        memoryCorrectionStandardError( ...
            :, iThreshold, iAmplitude) = ...
            result.memoryCorrectionStandardError;

        pCensored(iThreshold, iAmplitude) = ...
            result.pCensored;

        terminalProbability(iThreshold, iAmplitude) = ...
            result.F(end);

        terminalStandardError(iThreshold, iAmplitude) = ...
            result.terminalStandardError;
    end
end

family = struct();

family.x = envelope.x;
family.thresholdValues = thresholdValues;
family.amplitudeValues = amplitudeValues;

family.F = F;
family.S = S;
family.pLocal = pLocal;
family.memoryCorrection = memoryCorrection;

family.standardError = standardError;
family.memoryCorrectionStandardError = ...
    memoryCorrectionStandardError;

family.pCensored = pCensored;
family.terminalProbability = terminalProbability;
family.terminalStandardError = terminalStandardError;

end
