function [Areal, info, Acomplex] = buildAFromOWNS(solution, cfg)
%BUILDAFROMOWNS Construct real- and proper-complex energy matrices.
%
% Outputs:
%
%   Areal{n}
%       Real symmetric matrix for real Gaussian coefficients:
%
%           q = B*w,     w ~ N(0,I)
%           e = w.'*Areal*w
%
%       Areal = real(B'*W*B).
%
%   Acomplex{n}
%       Hermitian matrix for proper-complex Gaussian coefficients:
%
%           q = B*z,     z ~ CN(0,I)
%           e = z'*Acomplex*z
%
%       Acomplex = B'*W*B.
%
% MATLAB version: R2020b

Ny = numel(solution.eta);
Nx = size(solution.q, 2);
r = size(solution.q, 3);

numEnergyVariables = 5;

if isfield(cfg, 'owns') && ...
        isfield(cfg.owns, 'numEnergyVariables')
    numEnergyVariables = cfg.owns.numEnergyVariables;
end

energyDimension = numEnergyVariables * Ny;
energyIndices = 1:energyDimension;

if size(solution.gram_W, 1) ~= energyDimension
    error('buildAFromOWNS:WeightDimensionMismatch', ...
        ['gram_W has %d rows, while the configured energy state has ', ...
         '%d rows.'], ...
        size(solution.gram_W, 1), energyDimension);
end

Areal = cell(Nx, 1);
Acomplex = cell(Nx, 1);

symmetryDefectReal = zeros(Nx, 1);
hermitianDefectComplex = zeros(Nx, 1);

minimumEigenvalueReal = zeros(Nx, 1);
maximumEigenvalueReal = zeros(Nx, 1);

minimumEigenvalueComplex = zeros(Nx, 1);
maximumEigenvalueComplex = zeros(Nx, 1);

meanEnergy = zeros(Nx, 1);
varianceEnergyReal = zeros(Nx, 1);
varianceEnergyComplex = zeros(Nx, 1);

weightTolerance = 1e-12;

if isfield(cfg, 'owns') && ...
        isfield(cfg.owns, 'weightTolerance')
    weightTolerance = cfg.owns.weightTolerance;
end

fprintf('Constructing OWNS stochastic energy matrices:\n');

progressInterval = max(floor(Nx / 10), 1);

for n = 1:Nx
    Bn = squeeze( ...
        solution.q(energyIndices, n, :));

    weights = full(real( ...
        solution.gram_W(:, n)));

    smallNegative = ...
        weights < 0 & weights >= -weightTolerance;

    weights(smallNegative) = 0;

    if any(weights < 0)
        error('buildAFromOWNS:NegativeWeight', ...
            'A material negative energy weight occurs at station %d.', n);
    end

    % Energy-weighted Gram matrix.
    weightedB = bsxfun( ...
        @times, sqrt(weights), Bn);

    AcomplexN = weightedB' * weightedB;

    hermitianDefectComplex(n) = ...
        norm(AcomplexN - AcomplexN', 'fro') ...
        / max(norm(AcomplexN, 'fro'), eps);

    AcomplexN = ...
        0.5 * (AcomplexN + AcomplexN');

    % For real stochastic coefficients, the imaginary skew-symmetric
    % component contributes zero to w.'*A*w.
    ArealN = real(AcomplexN);

    symmetryDefectReal(n) = ...
        norm(ArealN - ArealN.', 'fro') ...
        / max(norm(ArealN, 'fro'), eps);

    ArealN = ...
        0.5 * (ArealN + ArealN.');

    lambdaReal = real(eig(ArealN));
    lambdaComplex = real(eig(AcomplexN));

    realScale = max(max(abs(lambdaReal)), eps);
    complexScale = max(max(abs(lambdaComplex)), eps);

    if min(lambdaReal) < ...
            -cfg.numerics.psdTol * realScale
        warning('buildAFromOWNS:RealPSDDefect', ...
            ['Real-coordinate A{%d} has minimum eigenvalue ', ...
             '%.6e.'], n, min(lambdaReal));
    end

    if min(lambdaComplex) < ...
            -cfg.numerics.psdTol * complexScale
        warning('buildAFromOWNS:ComplexPSDDefect', ...
            ['Proper-complex A{%d} has minimum eigenvalue ', ...
             '%.6e.'], n, min(lambdaComplex));
    end

    Areal{n} = ArealN;
    Acomplex{n} = AcomplexN;

    minimumEigenvalueReal(n) = min(lambdaReal);
    maximumEigenvalueReal(n) = max(lambdaReal);

    minimumEigenvalueComplex(n) = min(lambdaComplex);
    maximumEigenvalueComplex(n) = max(lambdaComplex);

    % Both conventions have the same mean:
    %
    %   E[e] = trace(A).
    meanEnergy(n) = real(trace(AcomplexN));

    % Real Gaussian:
    %
    %   Var(e) = 2*trace(Areal^2).
    varianceEnergyReal(n) = ...
        2 * trace(ArealN * ArealN);

    % Proper complex Gaussian:
    %
    %   Var(e) = trace(Acomplex^2).
    varianceEnergyComplex(n) = ...
        real(trace(AcomplexN * AcomplexN));

    if mod(n, progressInterval) == 0 || n == Nx
        fprintf('  completed station %d of %d\n', n, Nx);
    end
end

info = struct();

info.numStations = Nx;
info.stochasticRank = r;
info.energyDimension = energyDimension;
info.energyIndices = energyIndices;

info.meanEnergy = meanEnergy;
info.varianceEnergyReal = varianceEnergyReal;
info.varianceEnergyComplex = varianceEnergyComplex;

info.symmetryDefectReal = symmetryDefectReal;
info.hermitianDefectComplex = hermitianDefectComplex;

info.minimumEigenvalueReal = minimumEigenvalueReal;
info.maximumEigenvalueReal = maximumEigenvalueReal;

info.minimumEigenvalueComplex = minimumEigenvalueComplex;
info.maximumEigenvalueComplex = maximumEigenvalueComplex;

info.maximumSymmetryDefectReal = ...
    max(symmetryDefectReal);

info.maximumHermitianDefectComplex = ...
    max(hermitianDefectComplex);

end