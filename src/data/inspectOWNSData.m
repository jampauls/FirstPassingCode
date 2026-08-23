function report = inspectOWNSData(solution, cfg)
%INSPECTOWNSDATA Inspect and validate the OWNS solution structure.
%
%   report = inspectOWNSData(solution, cfg)
%
% This function does not modify or copy the large q array beyond ordinary
% MATLAB structure access.
%
% MATLAB version: R2020b

if nargin < 2
    cfg = struct();
end

requiredFields = { ...
    'q_in', ...
    'q', ...
    'gram_W', ...
    'x', ...
    'y', ...
    'eta'};

for j = 1:numel(requiredFields)
    if ~isfield(solution, requiredFields{j})
        error('inspectOWNSData:MissingField', ...
            'solution is missing required field "%s".', ...
            requiredFields{j});
    end
end

qSize = size(solution.q);
qInSize = size(solution.q_in);
weightSize = size(solution.gram_W);
xSize = size(solution.x);
ySize = size(solution.y);

Ny = numel(solution.eta);
Nx = qSize(2);
r = qSize(3);
nq = qSize(1);

expectedNq = 6 * Ny;
expectedEnergyDimension = 5 * Ny;

if nq ~= expectedNq
    error('inspectOWNSData:UnexpectedStateDimension', ...
        ['solution.q has %d rows, but 6*Ny = %d. ', ...
         'Check the state ordering.'], ...
        nq, expectedNq);
end

if qInSize(1) ~= expectedEnergyDimension
    error('inspectOWNSData:UnexpectedInletDimension', ...
        ['solution.q_in has %d rows, but 5*Ny = %d.'], ...
        qInSize(1), expectedEnergyDimension);
end

if qInSize(2) ~= r
    error('inspectOWNSData:InletRankMismatch', ...
        ['solution.q_in has %d columns, while solution.q has ', ...
         '%d stochastic columns.'], ...
        qInSize(2), r);
end

if any(weightSize ~= [expectedEnergyDimension, Nx])
    error('inspectOWNSData:WeightDimensionMismatch', ...
        ['solution.gram_W has size %d-by-%d; expected %d-by-%d.'], ...
        weightSize(1), weightSize(2), ...
        expectedEnergyDimension, Nx);
end

if size(solution.x, 2) ~= Nx || ...
        size(solution.y, 2) ~= Nx
    error('inspectOWNSData:GridStationMismatch', ...
        'The physical grid must contain Nx columns.');
end

if size(solution.x, 1) ~= Ny || ...
        size(solution.y, 1) ~= Ny
    error('inspectOWNSData:GridWallNormalMismatch', ...
        'The physical grid must contain Ny rows.');
end

% gram_W is small enough to inspect as a full vector, but avoid modifying
% the source structure itself.
weightValues = full(solution.gram_W(:));

weightMinimum = min(weightValues);
weightMaximum = max(weightValues);

clear weightValues;

weightTolerance = 1e-12;

if isfield(cfg, 'owns') && ...
        isfield(cfg.owns, 'weightTolerance')
    weightTolerance = cfg.owns.weightTolerance;
end

if weightMinimum < -weightTolerance
    error('inspectOWNSData:NegativeEnergyWeight', ...
        ['gram_W contains a negative entry %.6e that exceeds the ', ...
         'configured tolerance.'], weightMinimum);
end

% The first streamwise q slice is projected, so it need not equal q_in.
B0Projected = squeeze( ...
    solution.q(1:expectedEnergyDimension, 1, :));

projectionDifference = norm( ...
    B0Projected - solution.q_in, 'fro') ...
    / max(norm(solution.q_in, 'fro'), eps);

report = struct();

report.Nx = Nx;
report.Ny = Ny;
report.r = r;
report.nq = nq;
report.energyDimension = expectedEnergyDimension;

report.qSize = qSize;
report.qInSize = qInSize;
report.weightSize = weightSize;
report.xSize = xSize;
report.ySize = ySize;

report.qIsComplex = ~isreal(solution.q);
report.qInIsComplex = ~isreal(solution.q_in);
report.covarianceFactorOrientation = 'synthesis';
report.gaussianConvention = 'real';

report.minimumWeight = weightMinimum;
report.maximumWeight = weightMaximum;
report.projectedInletRelativeDifference = projectionDifference;

if isfield(solution, 'beta')
    report.beta = solution.beta;
end

if isfield(solution, 'w')
    report.omega = solution.w;
else
    error('inspectOWNSData:MissingFrequency', ...
        'The temporal frequency must be provided in solution.w.');
end

if isfield(solution, 'fp')
    report.forcingNorm = norm(solution.fp(:));
else
    report.forcingNorm = NaN;
end

verbose = true;

if isfield(cfg, 'owns') && ...
        isfield(cfg.owns, 'verbose')
    verbose = cfg.owns.verbose;
end

if verbose
    fprintf('\nOWNS data inspection:\n');
    fprintf('  Nx                         = %d\n', Nx);
    fprintf('  Ny                         = %d\n', Ny);
    fprintf('  physical state dimension  = %d\n', nq);
    fprintf('  energy-bearing dimension  = %d\n', ...
        expectedEnergyDimension);
    fprintf('  stochastic rank           = %d\n', r);
    fprintf('  q is complex              = %d\n', report.qIsComplex);
    fprintf('  q_in is complex           = %d\n', report.qInIsComplex);
    fprintf('  minimum energy weight     = %.6e\n', weightMinimum);
    fprintf('  maximum energy weight     = %.6e\n', weightMaximum);
    fprintf('  projected inlet difference = %.6e\n', ...
        projectionDifference);

    if isfield(report, 'omega')
        fprintf('  temporal frequency        = %.6e\n', report.omega);
    end

    if isfield(report, 'beta')
        fprintf('  spanwise wavenumber       = %.6e\n', report.beta);
    end

    if isfinite(report.forcingNorm)
        fprintf('  forcing norm              = %.6e\n', ...
            report.forcingNorm);
    end

    fprintf('\n');
end

end
