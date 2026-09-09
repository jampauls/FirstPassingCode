function precacheOWNSEnsemble(fileList, cfg)
%PRECACHEOWNSENSEMBLE Build missing or stale mode caches in parallel.

if isempty(fileList)
    return;
end

cacheDir = fullfile('data', 'reducedCache');
if isfield(cfg, 'owns') && isfield(cfg.owns, 'reducedCacheDirectory') && ...
        ~isempty(cfg.owns.reducedCacheDirectory)
    cacheDir = cfg.owns.reducedCacheDirectory;
end

if ~isfolder(cacheDir)
    mkdir(cacheDir);
end

numWorkers = localRecommendedWorkerCount(fileList, cfg);
useParallel = numWorkers > 1 && exist('parpool', 'file') == 2 && ...
    license('test', 'Distrib_Computing_Toolbox');

if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool('local', numWorkers);
        catch exception
            warning('precacheOWNSEnsemble:ParallelUnavailable', ...
                'Parallel cache build unavailable; using serial I/O. %s', ...
                exception.message);
            useParallel = false;
        end
    end
end

if useParallel
    fprintf('Pre-caching with %d workers.\n', numWorkers);
    parfor (i = 1:numel(fileList), numWorkers)
        loadOWNSModeCached(fileList{i}, cfg);
    end
else
    fprintf('Pre-caching serially.\n');
    for i = 1:numel(fileList)
        loadOWNSModeCached(fileList{i}, cfg);
    end
end

end

function numWorkers = localRecommendedWorkerCount(fileList, cfg)
%LOCALRECOMMENDEDWORKERCOUNT Limit raw-file loads by RAM and CPU capacity.

maxWorkers = inf;
if isfield(cfg, 'parallel') && isfield(cfg.parallel, 'maxWorkers') && ...
        ~isempty(cfg.parallel.maxWorkers)
    maxWorkers = cfg.parallel.maxWorkers;
end

sourceBytes = cellfun(@(file) dir(file), fileList, 'UniformOutput', false);
sourceBytes = cellfun(@(info) info.bytes, sourceBytes);
largestSourceBytes = max(sourceBytes);

numCores = feature('numcores');
cpuLimit = max(numCores - 2, 1);

memoryLimit = 1;
if isunix && isfile('/proc/meminfo')
    memInfo = fileread('/proc/meminfo');
    availableToken = regexp(memInfo, 'MemAvailable:\s+(\d+)', ...
        'tokens', 'once');
    totalToken = regexp(memInfo, 'MemTotal:\s+(\d+)', ...
        'tokens', 'once');
    availableKB = sscanf(availableToken{1}, '%f');
    totalKB = sscanf(totalToken{1}, '%f');

    reserveBytes = max(0.25 * totalKB * 1024, 4 * 1024^3);
    workerBytes = max(3 * largestSourceBytes, 2 * 1024^3);
    memoryLimit = max(floor((availableKB * 1024 - reserveBytes) / ...
        workerBytes), 1);
end

numWorkers = min([numel(fileList), cpuLimit, memoryLimit, maxWorkers]);
numWorkers = max(floor(numWorkers), 1);

end