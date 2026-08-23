function problem = loadOWNSProblem(source, cfg)
%LOADOWNSPROBLEM Load OWNS data and compress it to A(x).
%
%   problem = loadOWNSProblem(dataFile, cfg)
%   problem = loadOWNSProblem(solutionStruct, cfg)
%
% Output:
%   problem.x
%   problem.A
%   problem.eThresh
%   problem.meta
%
% The large physical factor q is not returned.
%
% MATLAB version: R2020b

if ischar(source) || isstring(source)
    dataFile = char(source);

    if ~isfile(dataFile)
        error('loadOWNSProblem:FileNotFound', ...
            'OWNS data file not found: %s', dataFile);
    end

    solutionVariable = 'solution';

    if isfield(cfg, 'owns') && ...
            isfield(cfg.owns, 'solutionVariable')
        solutionVariable = ...
            cfg.owns.solutionVariable;
    end

    fprintf('Loading OWNS solution from:\n  %s\n', dataFile);

    loaded = load(dataFile, solutionVariable);

    if ~isfield(loaded, solutionVariable)
        error('loadOWNSProblem:MissingSolutionVariable', ...
            'The file does not contain variable "%s".', ...
            solutionVariable);
    end

    solution = loaded.(solutionVariable);
    clear loaded;

    sourceDescription = dataFile;

elseif isstruct(source)
    solution = source;
    sourceDescription = 'in-memory solution structure';

else
    error('loadOWNSProblem:InvalidSource', ...
        'source must be a file path or solution structure.');
end

inspection = inspectOWNSData(solution, cfg);

[xCoordinate, coordinateMeta] = ...
    extractOWNSCoordinate( ...
        solution, cfg.coordinates);

[A, matrixInfo] = ...
    buildAFromOWNS(solution, cfg);

[eThresh, thresholdInfo] = ...
    buildTransitionThreshold( ...
        xCoordinate, A, cfg.threshold);

meta = struct();

meta.source = sourceDescription;
meta.dataType = 'OWNS';
meta.factorOrientation = 'synthesis';
meta.gaussianConvention = 'real';
meta.stateIsComplex = inspection.qIsComplex;

meta.Nx = inspection.Nx;
meta.Ny = inspection.Ny;
meta.nq = inspection.nq;
meta.r = inspection.r;
meta.energyDimension = inspection.energyDimension;

meta.stateOrdering = ...
    {'rho', 'u', 'v', 'w', 'T', 'p'};

meta.energyVariables = ...
    {'rho', 'u', 'v', 'w', 'T'};

meta.pressureIncludedInEnergy = false;
meta.initialProjectionApplied = true;
meta.distributedForcingIsZero = ...
    isfinite(inspection.forcingNorm) && ...
    inspection.forcingNorm == 0;

meta.inspection = inspection;
meta.coordinate = coordinateMeta;
meta.matrixInfo = matrixInfo;
meta.threshold = thresholdInfo;

if isfield(solution, 'beta')
    meta.beta = solution.beta;
end

if isfield(solution, 'w')
    meta.omega = solution.w;
else
    error('loadOWNSProblem:MissingFrequency', ...
        'The temporal frequency must be provided in solution.w.');
end

if isfield(solution, 'xi')
    meta.xi = solution.xi(:);
end

% Retain lightweight physical grid information for plotting.
meta.wallX = solution.x(1, :).';
meta.wallY = solution.y(1, :).';

if isfield(solution, 'z')
    meta.wallZ = solution.z(1, :).';
end

problem = struct();
problem.representation = 'energyMatrices';
problem.x = xCoordinate;
problem.A = A;
problem.eThresh = eThresh;
problem.meta = meta;

% The local solution variable and its large q field are released when this
% function returns.
clear solution;

end
