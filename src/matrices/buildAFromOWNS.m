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

[trimToDelta99, numDelta99, delta99, maskMethod, taperWidth] = ...
    getBoundaryLayerTrimOptions(solution, cfg, Nx);

if size(solution.gram_W, 1) ~= energyDimension
    error('buildAFromOWNS:WeightDimensionMismatch', ...
        ['gram_W has %d rows, while the configured energy state has ', ...
         '%d rows.'], ...
        size(solution.gram_W, 1), energyDimension);
end

if size(solution.gram_W, 1) ~= energyDimension
    error('buildAFromOWNS:WeightDimensionMismatch', ...
        ['gram_W has %d rows, while the configured energy state ', ...
         'has %d rows.'], ...
        size(solution.gram_W, 1), energyDimension);
end

if size(solution.gram_W, 2) ~= Nx
    error('buildAFromOWNS:WeightStationMismatch', ...
        ['gram_W has %d streamwise columns, while solution.q ', ...
         'contains %d stations.'], ...
        size(solution.gram_W, 2), Nx);
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

if trimToDelta99
    fprintf(['Restricting energy norm to %.6g delta99 ', ...
             'from the wall using %s masking.\n'], ...
        numDelta99, maskMethod);
end

progressInterval = max(floor(Nx / 10), 1);

for n = 1:Nx
    % ---------------------------------------------------------------------
    % Extract the propagated covariance synthesis factor in the physical
    % five-variable energy-bearing state:
    %
    %   [rho, u, v, w, T].
    %
    % solution.q also contains pressure as the sixth variable block.
    % Pressure is excluded from the energy norm by selecting only
    % energyIndices = 1:(5*Ny).
    %
    % Dimensions:
    %
    %   Bn : (5*Ny)-by-r
    % ---------------------------------------------------------------------

    Bn = squeeze( ...
        solution.q(energyIndices, n, :));

    % solution.gram_W is defined in the same physical variable space as
    % the selected rows of solution.q. It may be stored as sparse.
    %
    % Dimensions:
    %
    %   weights : (5*Ny)-by-1
    weights = full(real( ...
        solution.gram_W(:, n)));

    if trimToDelta99
        wallDistance = abs( ...
            solution.y(:, n) - solution.y(1, n));
        mask = computeDelta99Mask( ...
            wallDistance, numDelta99 * delta99(n), ...
            taperWidth * delta99(n), maskMethod);

        % The five energy variables occupy consecutive Ny-row blocks.
        weights = weights .* repmat(mask, numEnergyVariables, 1);
    end

    % ---------------------------------------------------------------------
    % Validate the diagonal energy weights
    % ---------------------------------------------------------------------

    smallNegative = ...
        weights < 0 ...
        & weights >= -weightTolerance;

    weights(smallNegative) = 0;

    if any(~isfinite(weights))
        error('buildAFromOWNS:NonfiniteWeight', ...
            ['A nonfinite energy weight occurs at streamwise ', ...
             'station %d.'], n);
    end

    if any(weights < 0)
        error('buildAFromOWNS:NegativeWeight', ...
            ['A material negative energy weight occurs at ', ...
             'streamwise station %d.'], n);
    end

    if all(weights == 0)
        error('buildAFromOWNS:ZeroEnergyMetric', ...
            ['All energy weights are zero at streamwise ', ...
             'station %d.'], n);
    end

    % ---------------------------------------------------------------------
    % Construct the energy-weighted physical factor
    %
    % Since W is diagonal,
    %
    %   W^(1/2) B
    %
    % is evaluated by rowwise multiplication.
    % ---------------------------------------------------------------------

    weightedB = bsxfun( ...
        @times, ...
        sqrt(weights), ...
        Bn);

    % ---------------------------------------------------------------------
    % Proper-complex stochastic-coordinate energy matrix
    %
    %   Acomplex = B^* W B
    %
    % This matrix is Hermitian positive semidefinite in exact arithmetic.
    % ---------------------------------------------------------------------

    AcomplexRaw = weightedB' * weightedB;

    complexScaleFro = max( ...
        norm(AcomplexRaw, 'fro'), eps);

    hermitianDefectComplex(n) = ...
        norm( ...
            AcomplexRaw - AcomplexRaw', ...
            'fro') ...
        / complexScaleFro;

    % Explicitly enforce Hermitian structure against roundoff.
    AcomplexN = ...
        0.5 * (AcomplexRaw + AcomplexRaw');

    % ---------------------------------------------------------------------
    % Real-Gaussian stochastic-coordinate energy matrix
    %
    % For real w,
    %
    %   e = w.' * Acomplex * w
    %     = w.' * real(Acomplex) * w,
    %
    % because the imaginary part of a Hermitian matrix is real
    % skew-symmetric and contributes zero to a real quadratic form.
    % ---------------------------------------------------------------------

    ArealRaw = real(AcomplexN);

    realScaleFro = max( ...
        norm(ArealRaw, 'fro'), eps);

    symmetryDefectReal(n) = ...
        norm( ...
            ArealRaw - ArealRaw.', ...
            'fro') ...
        / realScaleFro;

    % Explicitly enforce real symmetry against roundoff.
    ArealN = ...
        0.5 * (ArealRaw + ArealRaw.');

    % ---------------------------------------------------------------------
    % Positive-semidefinite diagnostics
    %
    % Do not clip matrix eigenvalues here. Report material negative
    % eigenvalues so that structural diagnostics are not hidden.
    % ---------------------------------------------------------------------

    lambdaReal = real(eig(ArealN));
    lambdaComplex = real(eig(AcomplexN));

    realScale = max( ...
        max(abs(lambdaReal)), eps);

    complexScale = max( ...
        max(abs(lambdaComplex)), eps);

    minimumEigenvalueReal(n) = ...
        min(lambdaReal);

    maximumEigenvalueReal(n) = ...
        max(lambdaReal);

    minimumEigenvalueComplex(n) = ...
        min(lambdaComplex);

    maximumEigenvalueComplex(n) = ...
        max(lambdaComplex);

    if minimumEigenvalueReal(n) < ...
            -cfg.numerics.psdTol * realScale

        warning('buildAFromOWNS:RealPSDDefect', ...
            ['Real-coordinate A{%d} has minimum eigenvalue ', ...
             '%.6e and relative PSD defect %.6e.'], ...
            n, ...
            minimumEigenvalueReal(n), ...
            max(-minimumEigenvalueReal(n), 0) / realScale);
    end

    if minimumEigenvalueComplex(n) < ...
            -cfg.numerics.psdTol * complexScale

        warning('buildAFromOWNS:ComplexPSDDefect', ...
            ['Proper-complex A{%d} has minimum eigenvalue ', ...
             '%.6e and relative PSD defect %.6e.'], ...
            n, ...
            minimumEigenvalueComplex(n), ...
            max(-minimumEigenvalueComplex(n), 0) ...
                / complexScale);
    end

    % ---------------------------------------------------------------------
    % Store the stochastic-space energy matrices
    % ---------------------------------------------------------------------

    Areal{n} = ArealN;
    Acomplex{n} = AcomplexN;

    % ---------------------------------------------------------------------
    % One-point energy moments
    %
    % Both conventions have the same mean:
    %
    %   E[e] = trace(A).
    %
    % Real Gaussian coefficients:
    %
    %   Var(e) = 2 trace(Areal^2).
    %
    % Proper-complex Gaussian coefficients:
    %
    %   Var(e) = trace(Acomplex^2).
    % ---------------------------------------------------------------------

    meanEnergy(n) = real(trace(AcomplexN));

    varianceEnergyReal(n) = real( ...
        2 * trace(ArealN * ArealN));

    varianceEnergyComplex(n) = real( ...
        trace(AcomplexN * AcomplexN));

    % ---------------------------------------------------------------------
    % Optional matrix consistency check
    %
    % The trace of Areal and Acomplex should agree because the imaginary
    % diagonal of a Hermitian matrix is zero.
    % ---------------------------------------------------------------------

    traceDifference = abs( ...
        trace(ArealN) ...
        - real(trace(AcomplexN)));

    traceScale = max( ...
        abs(real(trace(AcomplexN))), eps);

    if traceDifference / traceScale > 1e-12
        warning('buildAFromOWNS:TraceMismatch', ...
            ['Real and proper-complex energy matrices have a ', ...
             'relative trace mismatch %.6e at station %d.'], ...
            traceDifference / traceScale, n);
    end

    % ---------------------------------------------------------------------
    % Progress output
    % ---------------------------------------------------------------------

    if mod(n, progressInterval) == 0 || n == Nx
        fprintf('  completed station %d of %d\n', ...
            n, Nx);
    end
end

info = struct();

info.numStations = Nx;
info.stochasticRank = r;
info.energyDimension = energyDimension;
info.energyIndices = energyIndices;
info.boundaryLayerTrimEnabled = trimToDelta99;
info.numDelta99 = numDelta99;
info.boundaryLayerMaskMethod = maskMethod;
info.delta99TaperWidth = taperWidth;

if trimToDelta99
    info.delta99 = delta99;
end

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

% =========================================================================
function [enabled, numDelta99, delta99, maskMethod, taperWidth] = ...
    getBoundaryLayerTrimOptions(solution, cfg, Nx)

enabled = false;
numDelta99 = [];
delta99 = [];
maskMethod = 'binary';
taperWidth = 0;

if ~isfield(cfg, 'owns') || ...
        ~isfield(cfg.owns, 'trimToDelta99') || ...
        ~cfg.owns.trimToDelta99
    return;
end

enabled = true;

if ~isfield(cfg.owns, 'numDelta99') || ...
        isempty(cfg.owns.numDelta99)
    error('buildAFromOWNS:MissingNumDelta99', ...
        'cfg.owns.numDelta99 is required when trimToDelta99 is true.');
end

numDelta99 = cfg.owns.numDelta99;

validateattributes(numDelta99, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'cfg.owns.numDelta99');

if isfield(cfg.owns, 'delta99MaskMethod') && ...
        ~isempty(cfg.owns.delta99MaskMethod)
    maskMethod = lower(char(cfg.owns.delta99MaskMethod));
end

switch maskMethod
    case 'binary'
        taperWidth = 0;

    case 'cosinetaper'
        if isfield(cfg.owns, 'delta99TaperWidth') && ...
                ~isempty(cfg.owns.delta99TaperWidth)
            taperWidth = cfg.owns.delta99TaperWidth;
        else
            taperWidth = 0.25;
        end

        validateattributes(taperWidth, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'}, ...
            mfilename, 'cfg.owns.delta99TaperWidth');

    otherwise
        error('buildAFromOWNS:UnknownDelta99MaskMethod', ...
            ['cfg.owns.delta99MaskMethod must be ''binary'' or ', ...
             '''cosineTaper''.']);
end

if ~isfield(solution, 'delta') || isempty(solution.delta)
    error('buildAFromOWNS:MissingDelta99', ...
        ['solution.delta is required when trimToDelta99 is true. ', ...
         'It must contain one positive delta99 value per station.']);
end

delta99 = solution.delta(:);

if numel(delta99) ~= Nx
    error('buildAFromOWNS:Delta99LengthMismatch', ...
        ['solution.delta contains %d values; expected one delta99 ', ...
         'value for each of the %d streamwise stations.'], ...
        numel(delta99), Nx);
end

if any(~isfinite(delta99)) || any(delta99 <= 0)
    error('buildAFromOWNS:InvalidDelta99', ...
        'solution.delta must contain finite, strictly positive values.');
end

end

% =========================================================================
function mask = computeDelta99Mask( ...
    wallDistance, cutoffDistance, taperDistance, maskMethod)

switch maskMethod
    case 'binary'
        mask = double(wallDistance <= cutoffDistance);

    case 'cosinetaper'
        innerEdge = max(cutoffDistance - taperDistance, 0);
        outerEdge = cutoffDistance + taperDistance;

        mask = zeros(size(wallDistance));
        mask(wallDistance <= innerEdge) = 1;

        transition = wallDistance > innerEdge & ...
            wallDistance < outerEdge;

        phase = (wallDistance(transition) - innerEdge) ...
            / (outerEdge - innerEdge);

        mask(transition) = 0.5 * (1 + cos(pi * phase));

    otherwise
        error('buildAFromOWNS:UnknownDelta99MaskMethod', ...
            'Unsupported delta99 mask method "%s".', maskMethod);
end

end