function [Ufamily, diagnostics] = ...
    computeSequentialLoewnerEnvelope(A, options)
%COMPUTESEQUENTIALLOEWNERENVELOPE Build a deficit-based Loewner envelope.
%
% For an ordered symmetric matrix family A_k, construct
%
%   U_1 = A_1,
%
%   D_k = A_{k+1} - U_k,
%
%   U_{k+1} = U_k + (D_k)_+,
%
% where (D_k)_+ is the positive-semidefinite part of D_k.
%
% Then, up to floating-point error,
%
%   U_k >= A_j,  j <= k,
%
% and
%
%   U_{k+1} >= U_k.
%
% Options
% -------
% validationTolerance:
%   Relative tolerance used only for diagnostic pass/fail decisions.
%   Default: 1e-10.
%
% validateAllPrefixes:
%   If true, directly check U_k-A_j for all j<=k. Default: true.
%
% verbose:
%   Default: false.
%
% MATLAB version: R2020b

if nargin < 2
    options = struct();
end

if ~iscell(A) || isempty(A)
    error('computeSequentialLoewnerEnvelope:InvalidInput', ...
        'A must be a nonempty cell array.');
end

validationTolerance = getOption( ...
    options, 'validationTolerance', 1e-10);

validateAllPrefixes = getOption( ...
    options, 'validateAllPrefixes', true);

verbose = getOption(options, 'verbose', false);

Nx = numel(A);
r = size(A{1}, 1);

Ufamily = cell(Nx, 1);

minimumDeficitEigenvalue = nan(Nx - 1, 1);
maximumDeficitEigenvalue = nan(Nx - 1, 1);
positiveDeficitRank = zeros(Nx - 1, 1);

correctionTrace = zeros(Nx - 1, 1);
correctionFrobeniusNorm = zeros(Nx - 1, 1);
correctionSpectralNorm = zeros(Nx - 1, 1);
minimumCorrectionEigenvalue = zeros(Nx - 1, 1);

minimumCurrentDominanceEigenvalue = zeros(Nx, 1);
minimumPrefixDominanceEigenvalue = nan(Nx, 1);
minimumMonotoneIncrementEigenvalue = nan(Nx - 1, 1);

Ucurrent = symmetrizeMatrix(A{1}, r);
Ufamily{1} = Ucurrent;

minimumCurrentDominanceEigenvalue(1) = ...
    min(real(eig(Ucurrent - symmetrizeMatrix(A{1}, r))));

minimumPrefixDominanceEigenvalue(1) = ...
    minimumCurrentDominanceEigenvalue(1);

progressInterval = max(floor(Nx / 10), 1);

for k = 1:(Nx - 1)
    Anext = symmetrizeMatrix(A{k + 1}, r);

    deficit = Anext - Ucurrent;
    deficit = 0.5 * (deficit + deficit.');

    [V, Lambda] = eig(deficit);
    lambda = real(diag(Lambda));

    minimumDeficitEigenvalue(k) = min(lambda);
    maximumDeficitEigenvalue(k) = max(lambda);

    % Do not truncate positive eigenvalues using a modeling tolerance.
    % Retaining all floating-point-positive eigenvalues best preserves the
    % domination identity D_+ - D = (-D)_+.
    positiveLambda = max(lambda, 0);

    positiveDeficitRank(k) = nnz(lambda > 0);

    correction = ...
        V * diag(positiveLambda) * V.';

    correction = ...
        real(0.5 * (correction + correction.'));

    correctionTrace(k) = sum(positiveLambda);
    correctionFrobeniusNorm(k) = norm(correction, 'fro');
    correctionSpectralNorm(k) = max(positiveLambda);

    minimumCorrectionEigenvalue(k) = ...
        min(real(eig(correction)));

    Uprevious = Ucurrent;

    Ucurrent = Ucurrent + correction;
    Ucurrent = real(0.5 * (Ucurrent + Ucurrent.'));

    Ufamily{k + 1} = Ucurrent;

    currentDifference = Ucurrent - Anext;
    currentDifference = ...
        0.5 * (currentDifference + currentDifference.');

    minimumCurrentDominanceEigenvalue(k + 1) = ...
        min(real(eig(currentDifference)));

    monotoneDifference = Ucurrent - Uprevious;
    monotoneDifference = ...
        0.5 * (monotoneDifference + monotoneDifference.');

    minimumMonotoneIncrementEigenvalue(k) = ...
        min(real(eig(monotoneDifference)));

    if validateAllPrefixes
        minimumPrefixValue = Inf;

        for j = 1:(k + 1)
            Aj = symmetrizeMatrix(A{j}, r);

            difference = Ucurrent - Aj;
            difference = 0.5 * (difference + difference.');

            minimumPrefixValue = min( ...
                minimumPrefixValue, ...
                min(real(eig(difference))));
        end

        minimumPrefixDominanceEigenvalue(k + 1) = ...
            minimumPrefixValue;
    else
        minimumPrefixDominanceEigenvalue(k + 1) = ...
            minimumCurrentDominanceEigenvalue(k + 1);
    end

    if verbose && ...
            (mod(k + 1, progressInterval) == 0 || k + 1 == Nx)
        fprintf('  completed station %d of %d\n', k + 1, Nx);
    end
end

envelopeTrace = zeros(Nx, 1);
envelopeSpectralNorm = zeros(Nx, 1);
envelopeFrobeniusNorm = zeros(Nx, 1);

for k = 1:Nx
    Uk = Ufamily{k};

    envelopeTrace(k) = trace(Uk);
    envelopeSpectralNorm(k) = norm(Uk, 2);
    envelopeFrobeniusNorm(k) = norm(Uk, 'fro');
end

upperScale = max([ ...
    max(envelopeSpectralNorm), ...
    eps]);

minimumDominanceEigenvalue = ...
    min(minimumPrefixDominanceEigenvalue);

minimumMonotonicityEigenvalue = ...
    min(minimumMonotoneIncrementEigenvalue);

relativeDominanceDefect = ...
    max(-minimumDominanceEigenvalue, 0) ...
    / upperScale;

relativeMonotonicityDefect = ...
    max(-minimumMonotonicityEigenvalue, 0) ...
    / upperScale;

diagnostics = struct();

diagnostics.minimumDeficitEigenvalue = ...
    minimumDeficitEigenvalue;

diagnostics.maximumDeficitEigenvalue = ...
    maximumDeficitEigenvalue;

diagnostics.positiveDeficitRank = ...
    positiveDeficitRank;

diagnostics.correctionTrace = correctionTrace;
diagnostics.correctionFrobeniusNorm = ...
    correctionFrobeniusNorm;

diagnostics.correctionSpectralNorm = ...
    correctionSpectralNorm;

diagnostics.minimumCorrectionEigenvalue = ...
    minimumCorrectionEigenvalue;

diagnostics.minimumCurrentDominanceEigenvalue = ...
    minimumCurrentDominanceEigenvalue;

diagnostics.minimumPrefixDominanceEigenvalue = ...
    minimumPrefixDominanceEigenvalue;

diagnostics.minimumMonotoneIncrementEigenvalue = ...
    minimumMonotoneIncrementEigenvalue;

diagnostics.minimumDominanceEigenvalue = ...
    minimumDominanceEigenvalue;

diagnostics.minimumMonotonicityEigenvalue = ...
    minimumMonotonicityEigenvalue;

diagnostics.relativeDominanceDefect = ...
    relativeDominanceDefect;

diagnostics.relativeMonotonicityDefect = ...
    relativeMonotonicityDefect;

diagnostics.isNumericallyDominating = ...
    relativeDominanceDefect <= validationTolerance;

diagnostics.isNumericallyMonotone = ...
    relativeMonotonicityDefect <= validationTolerance;

diagnostics.envelopeTrace = envelopeTrace;
diagnostics.envelopeSpectralNorm = ...
    envelopeSpectralNorm;

diagnostics.envelopeFrobeniusNorm = ...
    envelopeFrobeniusNorm;

diagnostics.validateAllPrefixes = ...
    validateAllPrefixes;

diagnostics.validationTolerance = ...
    validationTolerance;

diagnostics.numStations = Nx;
diagnostics.stochasticDimension = r;

diagnostics.construction = ...
    'U_{k+1}=U_k+(A_{k+1}-U_k)_+';

end

% =========================================================================
function matrix = symmetrizeMatrix(matrix, expectedSize)

if ~isnumeric(matrix) || ...
        ~isequal(size(matrix), [expectedSize, expectedSize]) || ...
        any(~isfinite(matrix(:)))
    error('computeSequentialLoewnerEnvelope:InvalidMatrix', ...
        'Every member of A must be a finite square matrix.');
end

matrix = real(0.5 * (matrix + matrix.'));

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