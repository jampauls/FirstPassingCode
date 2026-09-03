function report = analyzeMatrixTrajectoryFeatures( ...
    matrixFamily, options)
%ANALYZEMATRIXTRAJECTORYFEATURES Analyze matrix-space compressibility.
%
% Let M_k be a family of real symmetric r-by-r matrices. Construct:
%
%   X(:,k) = svec(M_k).
%
% The singular values of X describe the intrinsic linear matrix-space
% dimension of the family.
%
% For candidate feature dimension m, the rank-m approximation is:
%
%   M_k^(m) = sum_{ell=1}^m theta(k,ell)*F_ell,
%
% where F_ell are Frobenius-orthonormal basis matrices.
%
% The function reports:
%
%   singular-value spectra;
%   cumulative Frobenius fractions;
%   numerical ranks;
%   coefficient trajectories;
%   Frobenius reconstruction errors;
%   spectral-norm reconstruction errors;
%   dimensions required for requested spectral error tolerances.
%
% Floating-point low rank is numerical evidence, not an exact algebraic
% certificate.
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
    error('analyzeMatrixTrajectoryFeatures:EmptyInput', ...
        'matrixFamily must contain at least one matrix.');
end

r = size(matrixFamily{1}, 1);

name = getOption(options, ...
    'name', 'matrix family');

candidateDimensions = getOption(options, ...
    'candidateDimensions', ...
    [1, 2, 3, 5, 10, 15, 20, 30, 40, 60, 80, 100]);

candidateDimensions = unique(round( ...
    candidateDimensions(:).'));

maximumPossibleRank = min( ...
    r * (r + 1) / 2, numMatrices);

candidateDimensions = candidateDimensions( ...
    candidateDimensions >= 0 ...
    & candidateDimensions <= maximumPossibleRank);

if isempty(candidateDimensions) || ...
        candidateDimensions(1) ~= 0
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

svecDimension = r * (r + 1) / 2;

% -------------------------------------------------------------------------
% Assemble the isometric matrix trajectory
% -------------------------------------------------------------------------

trajectory = zeros(svecDimension, numMatrices);

matrixFrobeniusNorm = zeros(numMatrices, 1);
matrixSpectralNorm = zeros(numMatrices, 1);
matrixTrace = zeros(numMatrices, 1);

minimumEigenvalue = zeros(numMatrices, 1);
maximumEigenvalue = zeros(numMatrices, 1);

for k = 1:numMatrices
    Mk = matrixFamily{k};

    if any(size(Mk) ~= [r, r])
        error('analyzeMatrixTrajectoryFeatures:DimensionMismatch', ...
            'matrixFamily{%d} has inconsistent dimensions.', k);
    end

    if ~isreal(Mk)
        error('analyzeMatrixTrajectoryFeatures:ComplexMatrix', ...
            'This diagnostic requires real symmetric matrices.');
    end

    Mk = 0.5 * (Mk + Mk.');

    trajectory(:, k) = ...
        symmetricMatrixToVector(Mk);

    matrixFrobeniusNorm(k) = norm(Mk, 'fro');
    matrixSpectralNorm(k) = norm(Mk, 2);
    matrixTrace(k) = trace(Mk);

    lambda = real(eig(Mk));

    minimumEigenvalue(k) = min(lambda);
    maximumEigenvalue(k) = max(lambda);
end

globalFrobeniusScale = max( ...
    matrixFrobeniusNorm);

globalSpectralScale = max( ...
    matrixSpectralNorm);

% -------------------------------------------------------------------------
% Compute the trajectory SVD through the smaller Gram matrix
%
%   Gram = trajectory' * trajectory.
%
% Its eigenvalues are squared singular values of trajectory.
% -------------------------------------------------------------------------

gramMatrix = trajectory.' * trajectory;
gramMatrix = 0.5 * (gramMatrix + gramMatrix.');

[Vright, Dgram] = eig(gramMatrix);

gramEigenvalues = real(diag(Dgram));

[gramEigenvalues, order] = ...
    sort(gramEigenvalues, 'descend');

Vright = real(Vright(:, order));

gramScale = max(max(abs(gramEigenvalues)), 1);

smallNegative = ...
    gramEigenvalues < 0 ...
    & gramEigenvalues >= -1e-12 * gramScale;

gramEigenvalues(smallNegative) = 0;

if any(gramEigenvalues < 0)
    warning('analyzeMatrixTrajectoryFeatures:GramPSDDefect', ...
        ['The trajectory Gram matrix has a material negative ', ...
         'eigenvalue %.6e.'], min(gramEigenvalues));

    gramEigenvalues = max(gramEigenvalues, 0);
end

singularValues = sqrt(gramEigenvalues);

singularValueEnergy = singularValues.^2;
totalSingularEnergy = sum(singularValueEnergy);

if totalSingularEnergy > 0
    cumulativeFrobeniusFraction = ...
        cumsum(singularValueEnergy) ...
        / totalSingularEnergy;
else
    cumulativeFrobeniusFraction = ...
        zeros(size(singularValueEnergy));
end

numericalRanks = zeros( ...
    numel(rankTolerances), 1);

for iTolerance = 1:numel(rankTolerances)
    numericalRanks(iTolerance) = ...
        relativeSingularRank( ...
            singularValues, ...
            rankTolerances(iTolerance));
end

% -------------------------------------------------------------------------
% Construct the nonzero left singular vectors:
%
%   U(:,ell) = X*V(:,ell)/sigma(ell).
% -------------------------------------------------------------------------

nonzeroTolerance = 1e-14;

if isempty(singularValues) || singularValues(1) == 0
    nonzeroRank = 0;
else
    nonzeroRank = nnz( ...
        singularValues ...
        > nonzeroTolerance * singularValues(1));
end

Uleft = zeros(svecDimension, nonzeroRank);

for ell = 1:nonzeroRank
    Uleft(:, ell) = ...
        trajectory * Vright(:, ell) ...
        / singularValues(ell);
end

% Coefficients:
%
%   coefficient(ell,k)
%       = <F_ell,M_k>_F
%       = sigma_ell * Vright(k,ell).
coefficientTrajectories = zeros( ...
    nonzeroRank, numMatrices);

if nonzeroRank > 0
    coefficientTrajectories = ...
        diag(singularValues(1:nonzeroRank)) ...
        * Vright(:, 1:nonzeroRank).';
end

% Store a limited number of basis matrices for later probability and
% convex-hull analyses.
basisCount = min( ...
    numBasisMatricesToStore, nonzeroRank);

basisMatrices = cell(basisCount, 1);

for ell = 1:basisCount
    basisMatrices{ell} = ...
        symmetricVectorToMatrix( ...
            Uleft(:, ell), r);
end

% -------------------------------------------------------------------------
% Candidate-dimension reconstruction diagnostics
%
% Residual is updated recursively:
%
%   R_m = R_{m-1} - u_m*(sigma_m*v_m').
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

maximumRelativeFrobeniusErrorLocal = zeros( ...
    numCandidates, 1);

medianRelativeFrobeniusErrorLocal = zeros( ...
    numCandidates, 1);

maximumRelativeFrobeniusErrorGlobal = zeros( ...
    numCandidates, 1);

maximumRelativeSpectralErrorLocal = zeros( ...
    numCandidates, 1);

medianRelativeSpectralErrorLocal = zeros( ...
    numCandidates, 1);

maximumRelativeSpectralErrorGlobal = zeros( ...
    numCandidates, 1);

residualTrajectory = trajectory;

currentFeatureDimension = 0;

for iCandidate = 1:numCandidates
    targetDimension = candidateDimensions(iCandidate);

    while currentFeatureDimension < targetDimension && ...
            currentFeatureDimension < nonzeroRank

        currentFeatureDimension = ...
            currentFeatureDimension + 1;

        residualTrajectory = ...
            residualTrajectory ...
            - Uleft(:, currentFeatureDimension) ...
            * coefficientTrajectories( ...
                currentFeatureDimension, :);
    end

    for k = 1:numMatrices
        residualVector = ...
            residualTrajectory(:, k);

        residualFrobeniusNorm = ...
            norm(residualVector);

        residualMatrix = ...
            symmetricVectorToMatrix( ...
                residualVector, r);

        residualSpectralNorm = ...
            norm(residualMatrix, 2);

        relativeFrobeniusErrorLocal(k, iCandidate) = ...
            residualFrobeniusNorm ...
            / max(matrixFrobeniusNorm(k), eps);

        relativeFrobeniusErrorGlobal(k, iCandidate) = ...
            residualFrobeniusNorm ...
            / max(globalFrobeniusScale, eps);

        relativeSpectralErrorLocal(k, iCandidate) = ...
            residualSpectralNorm ...
            / max(matrixSpectralNorm(k), eps);

        relativeSpectralErrorGlobal(k, iCandidate) = ...
            residualSpectralNorm ...
            / max(globalSpectralScale, eps);
    end

    maximumRelativeFrobeniusErrorLocal(iCandidate) = ...
        max(relativeFrobeniusErrorLocal(:, iCandidate));

    medianRelativeFrobeniusErrorLocal(iCandidate) = ...
        median(relativeFrobeniusErrorLocal(:, iCandidate));

    maximumRelativeFrobeniusErrorGlobal(iCandidate) = ...
        max(relativeFrobeniusErrorGlobal(:, iCandidate));

    maximumRelativeSpectralErrorLocal(iCandidate) = ...
        max(relativeSpectralErrorLocal(:, iCandidate));

    medianRelativeSpectralErrorLocal(iCandidate) = ...
        median(relativeSpectralErrorLocal(:, iCandidate));

    maximumRelativeSpectralErrorGlobal(iCandidate) = ...
        max(relativeSpectralErrorGlobal(:, iCandidate));
end

% -------------------------------------------------------------------------
% Feature dimensions required for target spectral errors
% -------------------------------------------------------------------------

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
% Package report
% -------------------------------------------------------------------------

report = struct();

report.name = name;

report.matrixDimension = r;
report.svecDimension = svecDimension;
report.numMatrices = numMatrices;

report.trajectory = trajectory;
report.gramMatrix = gramMatrix;

report.singularValues = singularValues;
report.singularValueEnergy = singularValueEnergy;
report.cumulativeFrobeniusFraction = ...
    cumulativeFrobeniusFraction;

report.rankTolerances = rankTolerances;
report.numericalRanks = numericalRanks;
report.nonzeroRank = nonzeroRank;

report.rightSingularVectors = Vright;
report.coefficientTrajectories = ...
    coefficientTrajectories;

report.numBasisMatricesStored = basisCount;
report.basisMatrices = basisMatrices;

report.matrixFrobeniusNorm = ...
    matrixFrobeniusNorm;

report.matrixSpectralNorm = ...
    matrixSpectralNorm;

report.matrixTrace = matrixTrace;
report.minimumEigenvalue = minimumEigenvalue;
report.maximumEigenvalue = maximumEigenvalue;

report.globalFrobeniusScale = ...
    globalFrobeniusScale;

report.globalSpectralScale = ...
    globalSpectralScale;

report.candidateDimensions = candidateDimensions;

report.relativeFrobeniusErrorLocal = ...
    relativeFrobeniusErrorLocal;

report.relativeFrobeniusErrorGlobal = ...
    relativeFrobeniusErrorGlobal;

report.relativeSpectralErrorLocal = ...
    relativeSpectralErrorLocal;

report.relativeSpectralErrorGlobal = ...
    relativeSpectralErrorGlobal;

report.maximumRelativeFrobeniusErrorLocal = ...
    maximumRelativeFrobeniusErrorLocal;

report.medianRelativeFrobeniusErrorLocal = ...
    medianRelativeFrobeniusErrorLocal;

report.maximumRelativeFrobeniusErrorGlobal = ...
    maximumRelativeFrobeniusErrorGlobal;

report.maximumRelativeSpectralErrorLocal = ...
    maximumRelativeSpectralErrorLocal;

report.medianRelativeSpectralErrorLocal = ...
    medianRelativeSpectralErrorLocal;

report.maximumRelativeSpectralErrorGlobal = ...
    maximumRelativeSpectralErrorGlobal;

report.spectralErrorTargets = ...
    spectralErrorTargets;

report.dimensionForGlobalSpectralTarget = ...
    dimensionForGlobalSpectralTarget;

report.dimensionForLocalSpectralTarget = ...
    dimensionForLocalSpectralTarget;

report.exactCertificateAvailable = false;
report.interpretation = ...
    ['Numerical linear matrix-space compression. Exactness requires ', ...
     'algebraically vanishing discarded singular values, which cannot ', ...
     'be certified from ordinary floating-point data alone.'];

if verbose
    fprintf('\nMatrix-trajectory feature analysis: %s\n', name);
    fprintf('  stochastic matrix size       = %d-by-%d\n', r, r);
    fprintf('  symmetric vector dimension   = %d\n', svecDimension);
    fprintf('  number of matrices           = %d\n', numMatrices);
    fprintf('  numerical nonzero rank       = %d\n', nonzeroRank);

    fprintf('\n  Numerical matrix-space ranks:\n');

    for iTolerance = 1:numel(rankTolerances)
        fprintf('    tolerance %.1e: %d\n', ...
            rankTolerances(iTolerance), ...
            numericalRanks(iTolerance));
    end

    fprintf('\n  Candidate reconstruction errors:\n');
    fprintf([ ...
        '    %7s %12s %12s %12s %12s\n'], ...
        'm', ...
        'Frob frac', ...
        'max spec/G', ...
        'max spec/L', ...
        'median spec/L');

    for iCandidate = 1:numCandidates
        m = candidateDimensions(iCandidate);

        if m == 0
            fraction = 0;
        elseif m <= numel(cumulativeFrobeniusFraction)
            fraction = cumulativeFrobeniusFraction(m);
        else
            fraction = 1;
        end

        fprintf([ ...
            '    %7d %12.5e %12.5e %12.5e %12.5e\n'], ...
            m, ...
            fraction, ...
            maximumRelativeSpectralErrorGlobal(iCandidate), ...
            maximumRelativeSpectralErrorLocal(iCandidate), ...
            medianRelativeSpectralErrorLocal(iCandidate));
    end

    fprintf('\n  Dimensions required for spectral error targets:\n');

    for iTarget = 1:numel(spectralErrorTargets)
        fprintf([ ...
            '    target %.1e: global=%s, local=%s\n'], ...
            spectralErrorTargets(iTarget), ...
            integerOrUnavailable( ...
                dimensionForGlobalSpectralTarget(iTarget)), ...
            integerOrUnavailable( ...
                dimensionForLocalSpectralTarget(iTarget)));
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
