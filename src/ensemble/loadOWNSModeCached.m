function mode = loadOWNSModeCached(dataFile, cfg)
%LOADOWNSMODECACHED Load a minimized OWNS mode, building the cache if needed.
%
%   mode = loadOWNSModeCached(dataFile, cfg)
%
% Building {x, A, meta} from a raw OWNS solution requires reading the
% full propagated factor solution.q, which can exceed tens of gigabytes
% per file for high-omega cases. This wrapper builds that much smaller
% representation once via loadOWNSModeRaw, saves it to a persistent
% cache directory, and reads the cached copy on every later call
% (including later runs of the workflow, not just later calls in the
% same run) so repeated passes over the ensemble never re-read the raw
% file.
%
% The cache is invalidated automatically if the raw source file's
% modification time or size changes, so an updated raw file is rebuilt
% rather than silently served from a stale cache. Any raw file without a
% cache entry yet (e.g. newly added to the ensemble directory) is built
% on this call without disturbing existing cache entries for other files.
%
% Optional cfg fields:
%   cfg.owns.reducedCacheDirectory   directory for cached minimal files,
%                                    created if it does not exist and
%                                    left in place after the run so later
%                                    runs reuse it. Default:
%                                    'data/reducedCache'.
%
% MATLAB version: R2020b

cacheDir = fullfile('data', 'reducedCache');

if isfield(cfg, 'owns') && ...
        isfield(cfg.owns, 'reducedCacheDirectory') && ...
        ~isempty(cfg.owns.reducedCacheDirectory)
    cacheDir = cfg.owns.reducedCacheDirectory;
end

if ~isfolder(cacheDir)
    mkdir(cacheDir);
end

cacheFile = fullfile(cacheDir, localCacheFileName(dataFile, cfg));

sourceInfo = dir(dataFile);
if isempty(sourceInfo)
    error('loadOWNSModeCached:SourceNotFound', ...
        'OWNS data file not found: %s', dataFile);
end

if isfile(cacheFile)
    loaded = load(cacheFile, 'mode', 'sourceDatenum', 'sourceBytes');

    isStale = ~isfield(loaded, 'sourceDatenum') || ...
        ~isfield(loaded, 'sourceBytes') || ...
        loaded.sourceDatenum ~= sourceInfo.datenum || ...
        loaded.sourceBytes ~= sourceInfo.bytes;

    if ~isStale
        mode = loaded.mode;
        return;
    end

    fprintf('Reduced cache is stale, rebuilding: %s\n', dataFile);
else
    fprintf('Building reduced cache (first access): %s\n', dataFile);
end

mode = loadOWNSModeRaw(dataFile, cfg);

sourceDatenum = sourceInfo.datenum;
sourceBytes = sourceInfo.bytes;

save(cacheFile, 'mode', 'sourceDatenum', 'sourceBytes', '-v7.3');

end

function name = localCacheFileName(dataFile, cfg)
%LOCALCACHEFILENAME Cache filename encoding the settings that affect it.

[~, baseName] = fileparts(dataFile);

numEnergyVariables = 5;
if isfield(cfg, 'owns') && isfield(cfg.owns, 'numEnergyVariables')
    numEnergyVariables = cfg.owns.numEnergyVariables;
end

coordMethod = 'wallArcLength';
referenceIndex = 1;
includeZ = true;
zeroOrigin = true;

if isfield(cfg, 'coordinates')
    if isfield(cfg.coordinates, 'method')
        coordMethod = cfg.coordinates.method;
    end
    if isfield(cfg.coordinates, 'referenceIndex')
        referenceIndex = cfg.coordinates.referenceIndex;
    end
    if isfield(cfg.coordinates, 'includeZ')
        includeZ = cfg.coordinates.includeZ;
    end
    if isfield(cfg.coordinates, 'zeroOrigin')
        zeroOrigin = cfg.coordinates.zeroOrigin;
    end
end

% Encodes every setting that changes the cached A(x) or x so a changed
% configuration cannot silently load a stale cache.
name = sprintf('%s_ev%d_%s_ref%d_z%d_zo%d.mat', ...
    baseName, numEnergyVariables, coordMethod, ...
    referenceIndex, includeZ, zeroOrigin);

end
