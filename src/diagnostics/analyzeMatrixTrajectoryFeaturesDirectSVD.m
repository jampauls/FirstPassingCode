function report = analyzeMatrixTrajectoryFeaturesDirectSVD( ...
    matrixFamily, options)
%ANALYZEMATRIXTRAJECTORYFEATURESDIRECTSVD Direct matrix-trajectory SVD.
%
% For a real symmetric matrix family M_k, form:
%
%   X(:,k) = svec(M_k),
%
% where svec preserves the Frobenius inner product. Then compute:
%
%   X = U*S*V'.
%
% The columns of U correspond to Frobenius-orthonormal matrix features:
%
%   F_l = smat(U(:,l)).
%
% The coefficient trajectory is:
%
%   theta(l,k) = S(l,l)*V(k,l).
%
% This direct SVD is more reliable for small singular values than an
% eigendecomposition of X'*X.
%
% Options:
%   name
%   candidateDimensions
%   rankTolerances
%   spectralErrorTargets
%   numBasisMatricesToStore
%   verbose
%
% MATLAB version: R2020b

if nargin < 2
    options = struct();
end

numMatrices = numel(matrixFamily);

if numMatrices == 0
    error('analyzeMatrixTrajectoryFeaturesDirectSVD:EmptyInput', ...
        'matrixFamily must contain at least one matrix.');
end

r = size(matrixFamily{1}, 1);
svecDimension = r * (r + 1) / 2;

name = getOption(options, ...
    'name', 'matrix family');

candidateDimensions = getOption(options, ...
    'candidateDimensions', ...
    [0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 60, 80, 100]);

candidateDimensions = unique(round( ...
    candidateDimensions(:).'));

maximumPossibleRank = min(svecDimension, numMatrices);

candidateDimensions = candidateDimensions( ...
    candidateDimensions >= 0 ...
    & candidateDimensions <= maximumPossibleRank);

if isempty(candidateDimensions) || candidateDimensions(1) ~= 0
    candidateDimensions = [0, candidateDimensions];
end

if candidateDimensions(end) ~= maximumPossibleRank
    candidateDimensions = [ ...
        candidateDimensions, maximumPossibleRank];
end

rankTolerances = getOption(options, ...
    'rankTolerances', ...
    [1e-6, 1e-8, 1e-10, 1e-12, 1e-14]);

rankTolerances = rankTolerances(:).';

spectralErrorTargets = getOption(options, ...
    'spectralErrorTargets', ...
    [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4, 1e-6]);

spectralErrorTargets = spectralErrorTargets(:).';

numBasisMatricesToStore = getOption(options, ...
    'numBasisMatricesToStore', ...
    min(30, maximumPossibleRank));

numBasisMatricesToStore = min( ...
    max(round(numBasisMatricesToStore), 0), ...
    maximumPossibleRank);

verbose = getOption(options, 'verbose', true);

% -------------------------------------------------------------------------
% Assemble trajectory
% -------------------------------------------------------------------------

trajectory = zeros(svecDimension, numMatrices);

matrixFrobeniusNorm = zeros(numMatrices, 1);
matrixSpectralNorm = zeros(numMatrices, 1);

for k = 1:numMatrices
    Mk = matrixFamily{k};

    if any(size(Mk) ~= [r, r])
        error( ...
            'analyzeMatrixTrajectoryFeaturesDirectSVD:DimensionMismatch', ...
            'matrixFamily{%d} has inconsistent dimensions.', k);
    end

    if ~isreal(Mk)
        error( ...
            'analyzeMatrixTrajectoryFeaturesDirectSVD:ComplexMatrix', ...
            'The matrix family must be real symmetric.');
    end

    Mk = 0.5 * (Mk + Mk.');

    trajectory(:, k) = ...
        symmetricMatrixToVector(Mk);

    matrixFrobeniusNorm(k) = norm(Mk, 'fro');
    matrixSpectralNorm(k) = norm(Mk, 2);
end

globalFrobeniusScale = max(matrixFrobeniusNorm);
globalSpectralScale = max(matrixSpectralNorm);

% -------------------------------------------------------------------------
% Direct economy SVD
% -------------------------------------------------------------------------

if verbose
    fprintf('\nDirect matrix-trajectory SVD: %s\n', name);
    fprintf('  trajectory size = %d-by-%d\n', ...
        size(trajectory, 1), size(trajectory, 2));
end

[Uleft, S, Vright] = svd(trajectory, 'econ');

singularValues = diag(S);
singularValueEnergy = singularValues.^2;

if sum(singularValueEnergy) > 0
    cumulativeFrobeniusFraction = ...
        cumsum(singularValueEnergy) ...
        / sum(singularValueEnergy);
else
    cumulativeFrobeniusFraction = ...
        zeros(size(singularValueEnergy));
end

numericalRanks = zeros(numel(rankTolerances), 1);

for iTolerance = 1:numel(rankTolerances)
    numericalRanks(iTolerance) = ...
        relativeSingularRank( ...
            singularValues, ...
            rankTolerances(iTolerance));
end

nonzeroRank = relativeSingularRank( ...
    singularValues, 1e-14);

coefficientTrajectories = ...
    S * Vright.';

basisCount = min( ...
    numBasisMatricesToStore, size(Uleft, 2));

basisMatrices = cell(basisCount, 1);

for ell = 1:basisCount
    basisMatrices{ell} = ...
        symmetricVectorToMatrix( ...
            Uleft(:, ell), r);
end

% -------------------------------------------------------------------------
% Reconstruction diagnostics
% -------------------------------------------------------------------------

numCandidates = numel(candidateDimensions);

relativeFrobeniusErrorLocal = zeros( ...
    numMatrices, numCandidates);

relativeFrobeniusErrorGlobal = zeros( ...
    numMatrices, numCandidates);

relativeSpectralErrorLocal = zeros( ...
    numMatrices, numCandidates);

relativeSpectralErrorGlobal = zeros( ...
    numMatrices, numCandidates);

maximumRelativeSpectralErrorLocal = zeros( ...
    numCandidates, 1);

medianRelativeSpectralErrorLocal = zeros( ...
    numCandidates, 1);

maximumRelativeSpectralErrorGlobal = zeros( ...
    numCandidates, 1);

maximumRelativeFrobeniusErrorGlobal = zeros( ...
    numCandidates, 1);

residualTrajectory = trajectory;
currentDimension = 0;

for iCandidate = 1:numCandidates
    targetDimension = candidateDimensions(iCandidate);

    while currentDimension < targetDimension && ...
            currentDimension < size(Uleft, 2)

        currentDimension = currentDimension + 1;

        residualTrajectory = ...
            residualTrajectory ...
            - Uleft(:, currentDimension) ...
            * coefficientTrajectories(currentDimension, :);
    end

    for k = 1:numMatrices
        residualVector = residualTrajectory(:, k);

        residualFrobenius = norm(residualVector);

        residualMatrix = ...
            symmetricVectorToMatrix( ...
                residualVector, r);

        residualSpectral = norm( ...
            residualMatrix, 2);

        relativeFrobeniusErrorLocal(k, iCandidate) = ...
            residualFrobenius ...
            / max(matrixFrobeniusNorm(k), eps);

        relativeFrobeniusErrorGlobal(k, iCandidate) = ...
            residualFrobenius ...
            / max(globalFrobeniusScale, eps);

        relativeSpectralErrorLocal(k, iCandidate) = ...
            residualSpectral ...
            / max(matrixSpectralNorm(k), eps);

        relativeSpectralErrorGlobal(k, iCandidate) = ...
            residualSpectral ...
            / max(globalSpectralScale, eps);
    end

    maximumRelativeSpectralErrorLocal(iCandidate) = ...
        max(relativeSpectralErrorLocal(:, iCandidate));

    medianRelativeSpectralErrorLocal(iCandidate) = ...
        median(relativeSpectralErrorLocal(:, iCandidate));

    maximumRelativeSpectralErrorGlobal(iCandidate) = ...
        max(relativeSpectralErrorGlobal(:, iCandidate));

    maximumRelativeFrobeniusErrorGlobal(iCandidate) = ...
        max(relativeFrobeniusErrorGlobal(:, iCandidate));
end

dimensionForGlobalSpectralTarget = nan( ...
    numel(spectralErrorTargets), 1);

dimensionForLocalSpectralTarget = nan( ...
    numel(spectralErrorTargets), 1);

for iTarget = 1:numel(spectralErrorTargets)
    target = spectralErrorTargets(iTarget);

    indexGlobal = find( ...
        maximumRelativeSpectralErrorGlobal <= target, ...
        1, 'first');

    if ~isempty(indexGlobal)
        dimensionForGlobalSpectralTarget(iTarget) = ...
            candidateDimensions(indexGlobal);
    end

    indexLocal = find( ...
        maximumRelativeSpectralErrorLocal <= target, ...
        1, 'first');

    if ~isempty(indexLocal)
        dimensionForLocalSpectralTarget(iTarget) = ...
            candidateDimensions(indexLocal);
    end
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

report = struct();

report.name = name;

report.matrixDimension = r;
report.svecDimension = svecDimension;
report.numMatrices = numMatrices;

report.trajectory = trajectory;

report.leftSingularVectors = Uleft;
report.singularValues = singularValues;
report.rightSingularVectors = Vright;

report.singularValueEnergy = singularValueEnergy;
report.cumulativeFrobeniusFraction = ...
    cumulativeFrobeniusFraction;

report.rankTolerances = rankTolerances;
report.numericalRanks = numericalRanks;
report.nonzeroRank = nonzeroRank;

report.coefficientTrajectories = ...
    coefficientTrajectories;

report.numBasisMatricesStored = basisCount;
report.basisMatrices = basisMatrices;

report.matrixFrobeniusNorm = matrixFrobeniusNorm;
report.matrixSpectralNorm = matrixSpectralNorm;

report.globalFrobeniusScale = globalFrobeniusScale;
report.globalSpectralScale = globalSpectralScale;

report.candidateDimensions = candidateDimensions;

report.relativeFrobeniusErrorLocal = ...
    relativeFrobeniusErrorLocal;

report.relativeFrobeniusErrorGlobal = ...
    relativeFrobeniusErrorGlobal;

report.relativeSpectralErrorLocal = ...
    relativeSpectralErrorLocal;

report.relativeSpectralErrorGlobal = ...
    relativeSpectralErrorGlobal;

report.maximumRelativeSpectralErrorLocal = ...
    maximumRelativeSpectralErrorLocal;

report.medianRelativeSpectralErrorLocal = ...
    medianRelativeSpectralErrorLocal;

report.maximumRelativeSpectralErrorGlobal = ...
    maximumRelativeSpectralErrorGlobal;

report.maximumRelativeFrobeniusErrorGlobal = ...
    maximumRelativeFrobeniusErrorGlobal;

report.spectralErrorTargets = spectralErrorTargets;

report.dimensionForGlobalSpectralTarget = ...
    dimensionForGlobalSpectralTarget;

report.dimensionForLocalSpectralTarget = ...
    dimensionForLocalSpectralTarget;

report.exactCertificateAvailable = false;

if verbose
    fprintf('  direct numerical ranks:\n');

    for iTolerance = 1:numel(rankTolerances)
        fprintf('    tolerance %.1e: %d\n', ...
            rankTolerances(iTolerance), ...
            numericalRanks(iTolerance));
    end

    fprintf('  candidate dimensions for global spectral targets:\n');

    for iTarget = 1:numel(spectralErrorTargets)
        fprintf('    target %.1e: %s\n', ...
            spectralErrorTargets(iTarget), ...
            integerOrUnavailable( ...
                dimensionForGlobalSpectralTarget(iTarget)));
    end

    fprintf('\n');
end

end

% =========================================================================
function numericalRank = relativeSingularRank( ...
    singularValues, tolerance)

if isempty(singularValues) || singularValues(1) == 0
    numericalRank = 0;
else
    numericalRank = nnz( ...
        singularValues > tolerance * singularValues(1));
end

end

% =========================================================================
function output = integerOrUnavailable(value)

if isnan(value)
    output = 'not reached';
else
    output = sprintf('%d', round(value));
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