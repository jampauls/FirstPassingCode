function report = validateOWNSQuadraticConstruction( ...
    solution, A, options)
%VALIDATEOWNSQUADRATICCONSTRUCTION Validate A_k from the OWNS factor.
%
% For selected streamwise stations, this function checks:
%
%   G_k = B_k' * W_k * B_k
%   A_k = real(G_k)
%
% and verifies:
%
%   (B_k*w)' * W_k * (B_k*w)
%       = w.' * A_k * w
%
% for deterministic and fixed pseudorandom real vectors w.
%
% This validates a floating-point algebraic identity. It does not estimate
% a probability.
%
% Inputs:
%   solution   OWNS solution structure
%   A          real symmetric stochastic energy matrices
%   options    diagnostic options
%
% MATLAB version: R2020b

if nargin < 3
    options = struct();
end

Ny = numel(solution.eta);
Nx = size(solution.q, 2);
r = size(solution.q, 3);

energyDimension = 5 * Ny;
energyIndices = 1:energyDimension;

if numel(A) ~= Nx
    error('validateOWNSQuadraticConstruction:LengthMismatch', ...
        'A must contain one matrix per OWNS station.');
end

numValidationStations = getOption( ...
    options, 'numValidationStations', 10);

numBasisVectors = getOption( ...
    options, 'numBasisValidationVectors', 5);

numRandomVectors = getOption( ...
    options, 'numRandomValidationVectors', 5);

randomSeed = getOption( ...
    options, 'validationRandomSeed', 1729);

verbose = getOption(options, 'verbose', true);

numValidationStations = min( ...
    max(round(numValidationStations), 2), Nx);

stationIndices = unique(round( ...
    linspace(1, Nx, numValidationStations))).';

numBasisVectors = min( ...
    max(round(numBasisVectors), 0), r);

numRandomVectors = max(round(numRandomVectors), 0);

% -------------------------------------------------------------------------
% Diagnostic real coefficient vectors
% -------------------------------------------------------------------------

Wtest = zeros(r, numBasisVectors + numRandomVectors);

for j = 1:numBasisVectors
    Wtest(j, j) = 1;
end

rng(randomSeed, 'twister');

if numRandomVectors > 0
    Wtest(:, numBasisVectors + (1:numRandomVectors)) = ...
        randn(r, numRandomVectors);

    % Normalize random vectors to avoid scale-dependent absolute errors.
    randomNorms = sqrt(sum( ...
        Wtest(:, numBasisVectors + (1:numRandomVectors)).^2, 1));

    Wtest(:, numBasisVectors + (1:numRandomVectors)) = ...
        bsxfun( ...
            @rdivide, ...
            Wtest(:, numBasisVectors + (1:numRandomVectors)), ...
            randomNorms);
end

numVectors = size(Wtest, 2);
numStations = numel(stationIndices);

hermitianResidual = zeros(numStations, 1);
symmetryResidualRaw = zeros(numStations, 1);
minimumEigenvalue = zeros(numStations, 1);

energyAbsoluteResidual = zeros( ...
    numStations, numVectors);

energyRelativeResidual = zeros( ...
    numStations, numVectors);

matrixConstructionResidual = zeros(numStations, 1);

for iStation = 1:numStations
    n = stationIndices(iStation);

    Bn = squeeze( ...
        solution.q(energyIndices, n, :));

    weights = full(real( ...
        solution.gram_W(:, n)));

    if any(weights < 0)
        error('validateOWNSQuadraticConstruction:NegativeWeight', ...
            'A negative energy weight occurs at station %d.', n);
    end

    weightedB = bsxfun(@times, sqrt(weights), Bn);

    Graw = weightedB' * weightedB;

    Gscale = max(norm(Graw, 2), 1);

    hermitianResidual(iStation) = ...
        norm(Graw - Graw', 2) / Gscale;

    Araw = real(Graw);

    Ascale = max(norm(Araw, 2), 1);

    symmetryResidualRaw(iStation) = ...
        norm(Araw - Araw.', 2) / Ascale;

    Asym = 0.5 * (Araw + Araw.');

    minimumEigenvalue(iStation) = ...
        min(real(eig(Asym)));

    matrixConstructionResidual(iStation) = ...
        norm(A{n} - Asym, 2) ...
        / max(norm(Asym, 2), eps);

    for j = 1:numVectors
        w = Wtest(:, j);

        qEnergy = Bn * w;

        physicalEnergy = real( ...
            qEnergy' * (weights .* qEnergy));

        stochasticEnergy = real( ...
            w.' * A{n} * w);

        residual = abs( ...
            physicalEnergy - stochasticEnergy);

        energyAbsoluteResidual(iStation, j) = ...
            residual;

        energyRelativeResidual(iStation, j) = ...
            residual ...
            / max([ ...
                abs(physicalEnergy), ...
                abs(stochasticEnergy), ...
                eps]);
    end
end

report = struct();

report.stationIndices = stationIndices;
report.testVectors = Wtest;

report.hermitianResidual = hermitianResidual;
report.symmetryResidualRaw = symmetryResidualRaw;
report.minimumEigenvalue = minimumEigenvalue;

report.energyAbsoluteResidual = ...
    energyAbsoluteResidual;

report.energyRelativeResidual = ...
    energyRelativeResidual;

report.matrixConstructionResidual = ...
    matrixConstructionResidual;

report.maximumHermitianResidual = ...
    max(hermitianResidual);

report.maximumSymmetryResidualRaw = ...
    max(symmetryResidualRaw);

report.minimumEigenvalueOverCheckedStations = ...
    min(minimumEigenvalue);

report.maximumEnergyAbsoluteResidual = ...
    max(energyAbsoluteResidual(:));

report.maximumEnergyRelativeResidual = ...
    max(energyRelativeResidual(:));

report.maximumMatrixConstructionResidual = ...
    max(matrixConstructionResidual);

report.isExactAlgebraicCertificate = false;
report.certificateStatus = ...
    'floating-point identity validation only';

if verbose
    fprintf('\nOWNS quadratic construction validation:\n');
    fprintf('  checked stations                 = %d\n', ...
        numStations);
    fprintf('  checked vectors per station      = %d\n', ...
        numVectors);
    fprintf('  max Hermitian residual           = %.6e\n', ...
        report.maximumHermitianResidual);
    fprintf('  max raw symmetry residual        = %.6e\n', ...
        report.maximumSymmetryResidualRaw);
    fprintf('  minimum checked eigenvalue       = %.6e\n', ...
        report.minimumEigenvalueOverCheckedStations);
    fprintf('  max A reconstruction residual    = %.6e\n', ...
        report.maximumMatrixConstructionResidual);
    fprintf('  max energy identity abs residual = %.6e\n', ...
        report.maximumEnergyAbsoluteResidual);
    fprintf('  max energy identity rel residual = %.6e\n', ...
        report.maximumEnergyRelativeResidual);
    fprintf('  certificate status               = %s\n\n', ...
        report.certificateStatus);
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