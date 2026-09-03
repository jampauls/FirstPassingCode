function study = computeCoupledMatrixFeatureProbabilityStudy( ...
    x, A, featureReport, featureDimensions, threshold, ...
    angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTECOUPLEDMATRIXFEATUREPROBABILITYSTUDY
% Compare matrix-feature truncations against the full first-passage model.
%
% For each feature dimension m:
%
%   A_k^(m)
%       = sum_{ell=1}^m theta_{ell,k} F_ell.
%
% The full and truncated models use identical RQMC directions in every
% replicate. This provides paired estimates of:
%
%   F_m(x) - F_full(x);
%   false-positive transition probability;
%   false-negative transition probability.
%
% A reconstructed matrix can be slightly indefinite because matrix-space
% SVD truncation does not preserve positive semidefiniteness. The radial
% first-passage formula remains meaningful for the homogeneous surrogate:
%
%   e_k^(m)(R*u) = R^2 u.' A_k^(m) u.
%
% Directions whose historical maximum is nonpositive have zero transition
% probability for a positive threshold.
%
% Inputs:
%   x                  Nx-by-1 streamwise coordinate
%   A                  full real symmetric matrix family
%   featureReport      direct-SVD feature report
%   featureDimensions  tested feature dimensions
%   threshold          positive scalar threshold
%   angularCfg         RQMC configuration
%   maxDetectionCfg    grid or Hermite maximum configuration
%   numericsCfg        numerical settings
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(A) ~= Nx
    error('computeCoupledMatrixFeatureProbabilityStudy:LengthMismatch', ...
        'A must contain one matrix per streamwise station.');
end

r = size(A{1}, 1);

validateattributes(threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'threshold');

featureDimensions = unique(round(featureDimensions(:).'));

if any(featureDimensions < 1)
    error('computeCoupledMatrixFeatureProbabilityStudy:InvalidDimension', ...
        'All feature dimensions must be positive integers.');
end

maximumFeatureDimension = max(featureDimensions);

if maximumFeatureDimension > featureReport.numBasisMatricesStored
    error('computeCoupledMatrixFeatureProbabilityStudy:BasisUnavailable', ...
        ['The largest requested feature dimension is %d, but only %d ', ...
         'basis matrices are stored in featureReport.'], ...
        maximumFeatureDimension, ...
        featureReport.numBasisMatricesStored);
end

if maximumFeatureDimension > ...
        size(featureReport.coefficientTrajectories, 1)

    error('computeCoupledMatrixFeatureProbabilityStudy:CoefficientUnavailable', ...
        'Insufficient feature coefficient trajectories are available.');
end

if size(featureReport.coefficientTrajectories, 2) ~= Nx
    error('computeCoupledMatrixFeatureProbabilityStudy:StationMismatch', ...
        'featureReport and A have inconsistent station counts.');
end

if ~strcmpi(angularCfg.method, 'rqmc')
    error('computeCoupledMatrixFeatureProbabilityStudy:AngularMethod', ...
        'angularCfg.method must be ''rqmc''.');
end

K = angularCfg.numDirections;
J = angularCfg.numReplicates;
numFeatures = numel(featureDimensions);

validateattributes(J, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, 'angularCfg.numReplicates');

% -------------------------------------------------------------------------
% Reconstruct all feature-truncated matrix families
% -------------------------------------------------------------------------

Afeature = cell(numFeatures, 1);
Gfeature = cell(numFeatures, 1);
GprimeFeature = cell(numFeatures, 1);

matrixDiagnostics = cell(numFeatures, 1);

for iFeature = 1:numFeatures
    m = featureDimensions(iFeature);

    Afeature{iFeature} = ...
        reconstructMatrixFamilyFromFeatures( ...
            featureReport, m);

    Gfeature{iFeature} = normalizeMatrixFamily( ...
        Afeature{iFeature}, threshold);

    if ismember(lower(maxDetectionCfg.method), ...
            {'hermite', 'adaptive'})

        derivativeCfg = struct();
        derivativeCfg.numerics = numericsCfg;

        [GprimeFeature{iFeature}, ~] = ...
            buildGprimeFamily( ...
                x, ...
                Gfeature{iFeature}, ...
                struct(), ...
                derivativeCfg, ...
                []);
    else
        GprimeFeature{iFeature} = [];
    end

    matrixDiagnostics{iFeature} = ...
        diagnoseReconstructedFamily( ...
            A, Afeature{iFeature});
end

% -------------------------------------------------------------------------
% Normalize the full matrix family
% -------------------------------------------------------------------------

Gfull = normalizeMatrixFamily(A, threshold);

if ismember(lower(maxDetectionCfg.method), ...
        {'hermite', 'adaptive'})

    derivativeCfg = struct();
    derivativeCfg.numerics = numericsCfg;

    [GprimeFull, ~] = buildGprimeFamily( ...
        x, Gfull, struct(), derivativeCfg, []);
else
    GprimeFull = [];
end

% -------------------------------------------------------------------------
% Replicate storage
% -------------------------------------------------------------------------

FfullReplicate = zeros(Nx, J);
FlocalFullReplicate = zeros(Nx, J);

FfeatureReplicate = zeros(Nx, numFeatures, J);
FlocalFeatureReplicate = zeros(Nx, numFeatures, J);

falsePositiveReplicate = zeros(Nx, numFeatures, J);
falseNegativeReplicate = zeros(Nx, numFeatures, J);

maximumEnvelopeAbsoluteErrorReplicate = ...
    zeros(Nx, numFeatures, J);

weightedEnvelopeAbsoluteErrorReplicate = ...
    zeros(Nx, numFeatures, J);

cfgLocal = struct();
cfgLocal.numerics = numericsCfg;
cfgLocal.maxDetection = maxDetectionCfg;
cfgLocal.gaussianConvention = 'real';

% Preserve negative gains in truncated matrix models. Do not flood the
% output with warnings from small indefinite components.
cfgLocal.numerics.clipSmallNegativeGains = false;

fprintf('\nCoupled matrix-feature probability study:\n');
fprintf('  full stochastic dimension = %d\n', r);
fprintf('  directions/replicate      = %d\n', K);
fprintf('  replicates                = %d\n', J);
fprintf('  feature dimensions        = %s\n\n', ...
    mat2str(featureDimensions));

for jReplicate = 1:J
    fprintf('  replicate %d of %d\n', ...
        jReplicate, J);

    replicateCfg = angularCfg;
    replicateCfg.replicateIndex = jReplicate;

    [U, weights] = generateAngularRule( ...
        r, replicateCfg);

    validateAngularRule(U, weights);

    % ---------------------------------------------------------------------
    % Full model
    % ---------------------------------------------------------------------

    fullMaxData = evaluateMaximumData( ...
        x, ...
        Gfull, ...
        GprimeFull, ...
        U, ...
        cfgLocal);

    fullProbability = computeTransitionCDF( ...
        fullMaxData.m, weights, r, cfgLocal);

    fullLocal = computeLocalTailAllowIndefinite( ...
        fullMaxData.a, weights, r);

    FfullReplicate(:, jReplicate) = ...
        fullProbability.F;

    FlocalFullReplicate(:, jReplicate) = ...
        fullLocal;

    % ---------------------------------------------------------------------
    % Feature models
    % ---------------------------------------------------------------------

    for iFeature = 1:numFeatures
        featureMaxData = evaluateMaximumData( ...
            x, ...
            Gfeature{iFeature}, ...
            GprimeFeature{iFeature}, ...
            U, ...
            cfgLocal);

        featureProbability = computeTransitionCDF( ...
            featureMaxData.m, weights, r, cfgLocal);

        featureLocal = computeLocalTailAllowIndefinite( ...
            featureMaxData.a, weights, r);

        FfeatureReplicate(:, iFeature, jReplicate) = ...
            featureProbability.F;

        FlocalFeatureReplicate(:, iFeature, jReplicate) = ...
            featureLocal;

        % The directional transition-tail values use the same directions
        % and therefore permit a paired false-positive/negative analysis.
        fullDirectionalTail = directionalFirstPassageTail( ...
            fullMaxData.m, r);

        featureDirectionalTail = directionalFirstPassageTail( ...
            featureMaxData.m, r);

        directionalDifference = ...
            featureDirectionalTail ...
            - fullDirectionalTail;

        falsePositiveReplicate(:, iFeature, jReplicate) = ...
            weights.' * max(directionalDifference, 0);

        falseNegativeReplicate(:, iFeature, jReplicate) = ...
            weights.' * max(-directionalDifference, 0);

        envelopeAbsoluteError = abs( ...
            featureMaxData.m - fullMaxData.m);

        maximumEnvelopeAbsoluteErrorReplicate( ...
            :, iFeature, jReplicate) = ...
            max(envelopeAbsoluteError, [], 1).';

        weightedEnvelopeAbsoluteErrorReplicate( ...
            :, iFeature, jReplicate) = ...
            weights.' * envelopeAbsoluteError;
    end
end

% -------------------------------------------------------------------------
% Replicate means and paired uncertainties
% -------------------------------------------------------------------------

Ffull = mean(FfullReplicate, 2);
FfullStandardError = ...
    std(FfullReplicate, 0, 2) / sqrt(J);

FlocalFull = mean(FlocalFullReplicate, 2);
FlocalFullStandardError = ...
    std(FlocalFullReplicate, 0, 2) / sqrt(J);

Ffeature = mean(FfeatureReplicate, 3);
FlocalFeature = mean(FlocalFeatureReplicate, 3);

FfeatureStandardError = ...
    std(FfeatureReplicate, 0, 3) / sqrt(J);

FlocalFeatureStandardError = ...
    std(FlocalFeatureReplicate, 0, 3) / sqrt(J);

cdfDifferenceReplicate = ...
    FfeatureReplicate ...
    - reshape(FfullReplicate, Nx, 1, J);

localDifferenceReplicate = ...
    FlocalFeatureReplicate ...
    - reshape(FlocalFullReplicate, Nx, 1, J);

cdfDifference = mean( ...
    cdfDifferenceReplicate, 3);

localDifference = mean( ...
    localDifferenceReplicate, 3);

cdfDifferenceStandardError = ...
    std(cdfDifferenceReplicate, 0, 3) / sqrt(J);

localDifferenceStandardError = ...
    std(localDifferenceReplicate, 0, 3) / sqrt(J);

falsePositiveProbability = ...
    mean(falsePositiveReplicate, 3);

falseNegativeProbability = ...
    mean(falseNegativeReplicate, 3);

falsePositiveStandardError = ...
    std(falsePositiveReplicate, 0, 3) / sqrt(J);

falseNegativeStandardError = ...
    std(falseNegativeReplicate, 0, 3) / sqrt(J);

maximumEnvelopeAbsoluteError = ...
    mean(maximumEnvelopeAbsoluteErrorReplicate, 3);

weightedEnvelopeAbsoluteError = ...
    mean(weightedEnvelopeAbsoluteErrorReplicate, 3);

% -------------------------------------------------------------------------
% Complete-CDF summary values
% -------------------------------------------------------------------------

maximumAbsoluteCdfDifference = zeros(numFeatures, 1);
maximumAbsoluteCdfDifferenceStandardError = zeros(numFeatures, 1);
maximumAbsoluteCdfDifferenceIndex = zeros(numFeatures, 1);

maximumAbsoluteLocalDifference = zeros(numFeatures, 1);
maximumFalsePositiveProbability = zeros(numFeatures, 1);
maximumFalseNegativeProbability = zeros(numFeatures, 1);

terminalCdfDifference = zeros(numFeatures, 1);
terminalCdfDifferenceStandardError = zeros(numFeatures, 1);

for iFeature = 1:numFeatures
    [maximumAbsoluteCdfDifference(iFeature), index] = ...
        max(abs(cdfDifference(:, iFeature)));

    maximumAbsoluteCdfDifferenceIndex(iFeature) = index;

    maximumAbsoluteCdfDifferenceStandardError(iFeature) = ...
        cdfDifferenceStandardError(index, iFeature);

    maximumAbsoluteLocalDifference(iFeature) = ...
        max(abs(localDifference(:, iFeature)));

    maximumFalsePositiveProbability(iFeature) = ...
        max(falsePositiveProbability(:, iFeature));

    maximumFalseNegativeProbability(iFeature) = ...
        max(falseNegativeProbability(:, iFeature));

    terminalCdfDifference(iFeature) = ...
        cdfDifference(end, iFeature);

    terminalCdfDifferenceStandardError(iFeature) = ...
        cdfDifferenceStandardError(end, iFeature);
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

study = struct();

study.x = x;
study.threshold = threshold;
study.fullStochasticDimension = r;

study.featureDimensions = featureDimensions;

study.FfullReplicate = FfullReplicate;
study.Ffull = Ffull;
study.FfullStandardError = FfullStandardError;

study.FlocalFullReplicate = FlocalFullReplicate;
study.FlocalFull = FlocalFull;
study.FlocalFullStandardError = ...
    FlocalFullStandardError;

study.FfeatureReplicate = FfeatureReplicate;
study.Ffeature = Ffeature;
study.FfeatureStandardError = ...
    FfeatureStandardError;

study.FlocalFeatureReplicate = ...
    FlocalFeatureReplicate;

study.FlocalFeature = FlocalFeature;
study.FlocalFeatureStandardError = ...
    FlocalFeatureStandardError;

study.cdfDifferenceReplicate = ...
    cdfDifferenceReplicate;

study.cdfDifference = cdfDifference;
study.cdfDifferenceStandardError = ...
    cdfDifferenceStandardError;

study.localDifferenceReplicate = ...
    localDifferenceReplicate;

study.localDifference = localDifference;
study.localDifferenceStandardError = ...
    localDifferenceStandardError;

study.falsePositiveProbability = ...
    falsePositiveProbability;

study.falseNegativeProbability = ...
    falseNegativeProbability;

study.falsePositiveStandardError = ...
    falsePositiveStandardError;

study.falseNegativeStandardError = ...
    falseNegativeStandardError;

study.maximumEnvelopeAbsoluteError = ...
    maximumEnvelopeAbsoluteError;

study.weightedEnvelopeAbsoluteError = ...
    weightedEnvelopeAbsoluteError;

study.maximumAbsoluteCdfDifference = ...
    maximumAbsoluteCdfDifference;

study.maximumAbsoluteCdfDifferenceStandardError = ...
    maximumAbsoluteCdfDifferenceStandardError;

study.maximumAbsoluteCdfDifferenceIndex = ...
    maximumAbsoluteCdfDifferenceIndex;

study.maximumAbsoluteLocalDifference = ...
    maximumAbsoluteLocalDifference;

study.maximumFalsePositiveProbability = ...
    maximumFalsePositiveProbability;

study.maximumFalseNegativeProbability = ...
    maximumFalseNegativeProbability;

study.terminalCdfDifference = ...
    terminalCdfDifference;

study.terminalCdfDifferenceStandardError = ...
    terminalCdfDifferenceStandardError;

study.matrixDiagnostics = matrixDiagnostics;
study.Afeature = Afeature;

study.numDirections = K;
study.numReplicates = J;
study.maximumMethod = maxDetectionCfg.method;

end

% =========================================================================
function G = normalizeMatrixFamily(A, threshold)

Nx = numel(A);
G = cell(Nx, 1);

for n = 1:Nx
    Gn = A{n} / threshold;
    G{n} = 0.5 * (Gn + Gn.');
end

end

% =========================================================================
function maxData = evaluateMaximumData( ...
    x, G, Gprime, U, cfgLocal)

switch lower(cfgLocal.maxDetection.method)

    case 'grid'
        maxData = computeRunningMaxGrid( ...
            G, U, cfgLocal);

    case 'hermite'
        maxData = computeRunningMaxHermite( ...
            x, G, Gprime, U, cfgLocal);

    case 'adaptive'
        maxData = computeRunningMaxAdaptive( ...
            x, G, Gprime, U, cfgLocal);

    otherwise
        error( ...
            'computeCoupledMatrixFeatureProbabilityStudy:MaximumMethod', ...
            'Unknown maximum method "%s".', ...
            cfgLocal.maxDetection.method);
end

end

% =========================================================================
function tail = directionalFirstPassageTail( ...
    runningMaximum, stochasticDimension)

tail = zeros(size(runningMaximum));

positive = runningMaximum > 0;

tail(positive) = gammainc( ...
    1 ./ (2 * runningMaximum(positive)), ...
    stochasticDimension / 2, ...
    'upper');

end

% =========================================================================
function pLocal = computeLocalTailAllowIndefinite( ...
    gain, weights, stochasticDimension)

[K, Nx] = size(gain);
weights = weights(:);

if numel(weights) ~= K
    error('computeLocalTailAllowIndefinite:WeightMismatch', ...
        'Weights are incompatible with the directional gains.');
end

pLocal = zeros(Nx, 1);

for n = 1:Nx
    gn = gain(:, n);

    positive = gn > 0;

    directionalTail = zeros(K, 1);

    directionalTail(positive) = gammainc( ...
        1 ./ (2 * gn(positive)), ...
        stochasticDimension / 2, ...
        'upper');

    pLocal(n) = weights.' * directionalTail;
end

end

% =========================================================================
function diagnostics = diagnoseReconstructedFamily( ...
    Afull, Aapprox)

Nx = numel(Afull);

minimumEigenvalue = zeros(Nx, 1);
relativePsdDefect = zeros(Nx, 1);

relativeSpectralErrorLocal = zeros(Nx, 1);
relativeSpectralErrorGlobal = zeros(Nx, 1);

globalScale = max(cellfun( ...
    @(M) norm(M, 2), Afull));

for n = 1:Nx
    Ahat = 0.5 * (Aapprox{n} + Aapprox{n}.');
    Atrue = 0.5 * (Afull{n} + Afull{n}.');

    lambda = real(eig(Ahat));

    minimumEigenvalue(n) = min(lambda);

    scale = max(max(abs(lambda)), eps);

    relativePsdDefect(n) = ...
        max(-minimumEigenvalue(n), 0) / scale;

    residualNorm = norm(Ahat - Atrue, 2);

    relativeSpectralErrorLocal(n) = ...
        residualNorm / max(norm(Atrue, 2), eps);

    relativeSpectralErrorGlobal(n) = ...
        residualNorm / max(globalScale, eps);
end

diagnostics = struct();

diagnostics.minimumEigenvalue = minimumEigenvalue;
diagnostics.relativePsdDefect = relativePsdDefect;

diagnostics.relativeSpectralErrorLocal = ...
    relativeSpectralErrorLocal;

diagnostics.relativeSpectralErrorGlobal = ...
    relativeSpectralErrorGlobal;

diagnostics.minimumEigenvalueOverX = ...
    min(minimumEigenvalue);

diagnostics.maximumRelativePsdDefect = ...
    max(relativePsdDefect);

diagnostics.numIndefiniteStations = ...
    nnz(relativePsdDefect > 1e-10);

diagnostics.maximumRelativeSpectralErrorLocal = ...
    max(relativeSpectralErrorLocal);

diagnostics.maximumRelativeSpectralErrorGlobal = ...
    max(relativeSpectralErrorGlobal);

end
