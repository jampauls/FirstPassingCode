function mode = loadOWNSModeRaw(dataFile, cfg)
%LOADOWNSMODERAW Load one OWNS (omega,beta) mode without a threshold.
%
%   mode = loadOWNSModeRaw(dataFile, cfg)
%
% This mirrors loadOWNSProblem, but stops short of building a
% single-mode transition threshold. It is intended for ensemble
% processing, where the transition threshold is defined from the total
% energy summed across every mode.
%
% The raw streamwise grid is downsampled (see downsampleOWNSStations) via
% cfg.owns.numCacheStations/cfg.owns.cacheXCutoffCoefficient before any
% per-station energy matrix is built, so the expensive part of this
% function only runs at the reduced set of stations that are actually
% cached.
%
% Output fields:
%   mode.x      streamwise coordinate for this mode (possibly downsampled)
%   mode.A      cell array of real energy matrices A_i(x)
%   mode.meta   metadata, including meta.omega, meta.beta,
%               meta.numRawStations, and meta.numCacheStations
%
% MATLAB version: R2020b

if ~isfile(dataFile)
    error('loadOWNSModeRaw:FileNotFound', ...
        'OWNS data file not found: %s', dataFile);
end

solutionVariable = 'solution';

if isfield(cfg, 'owns') && isfield(cfg.owns, 'solutionVariable')
    solutionVariable = cfg.owns.solutionVariable;
end

loaded = load(dataFile, solutionVariable);

if ~isfield(loaded, solutionVariable)
    error('loadOWNSModeRaw:MissingSolutionVariable', ...
        'The file does not contain variable "%s": %s', ...
        solutionVariable, dataFile);
end

solution = loaded.(solutionVariable);
clear loaded;

numRawStations = size(solution.q, 2);

solution = downsampleOWNSStations(solution, cfg);

inspection = inspectOWNSData(solution, cfg);

[xCoordinate, coordinateMeta] = ...
    extractOWNSCoordinate(solution, cfg.coordinates);

[A, matrixInfo] = buildAFromOWNS(solution, cfg);

meta = struct();
meta.source = dataFile;
meta.nq = inspection.nq;
meta.r = inspection.r;
meta.Ny = inspection.Ny;
meta.energyDimension = inspection.energyDimension;
meta.coordinate = coordinateMeta;
meta.matrixInfo = matrixInfo;
meta.numRawStations = numRawStations;
meta.numCacheStations = inspection.Nx;

if isfield(solution, 'w')
    meta.omega = solution.w;
else
    error('loadOWNSModeRaw:MissingFrequency', ...
        'solution.w (temporal frequency) is required: %s', dataFile);
end

if isfield(solution, 'beta')
    meta.beta = solution.beta;
else
    meta.beta = NaN;
end

mode = struct();
mode.x = xCoordinate;
mode.A = A;
mode.meta = meta;

clear solution;

end
