% runOWNSCacheStalenessTest.m
%
% Verifies the reduced-cache logic in loadOWNSModeCached:
%
%   1. the cache directory is persistent: an unchanged raw file is never
%      rebuilt on a later call;
%   2. a raw file whose modification time/size changes is detected as
%      stale and rebuilt, without touching other cache entries;
%   3. a raw file with no cache entry yet (e.g. newly added to the
%      ensemble directory) is built on demand, without rebuilding
%      already-cached files.
%
% Uses copies of two existing OWNS data files so the persistent
% data/reducedCache directory and the original data files are untouched.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS reduced-cache staleness test\n');
fprintf('============================================================\n\n');

sourceDataDir = fullfile(projectRoot, 'data');
tempDataDir = fullfile(tempdir, 'FirstPassingCodeCacheStalenessData');
cacheDir = fullfile(tempdir, 'FirstPassingCodeCacheStalenessCache');

if isfolder(tempDataDir)
    rmdir(tempDataDir, 's');
end
if isfolder(cacheDir)
    rmdir(cacheDir, 's');
end
mkdir(tempDataDir);

sourceFiles = discoverOWNSEnsembleFiles(sourceDataDir);
if numel(sourceFiles) < 2
    error('runOWNSCacheStalenessTest:TooFewFiles', ...
        'Expected at least 2 OWNS mode files in data/, found %d.', ...
        numel(sourceFiles));
end

[~, name1, ext1] = fileparts(sourceFiles{1});
[~, name2, ext2] = fileparts(sourceFiles{2});

file1 = fullfile(tempDataDir, [name1, ext1]);
file2 = fullfile(tempDataDir, [name2, ext2]);

copyfile(sourceFiles{1}, file1);
copyfile(sourceFiles{2}, file2);

%% Configuration (only owns/coordinates fields matter for this test)

cfg = struct();
cfg.owns = struct();
cfg.owns.reducedCacheDirectory = cacheDir;
cfg.owns.solutionVariable = 'solution';
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;
cfg.owns.validateWeights = true;
cfg.owns.weightTolerance = 1e-12;
cfg.owns.verbose = false;

cfg.coordinates = struct();
cfg.coordinates.method = 'wallArcLength';
cfg.coordinates.referenceIndex = 1;
cfg.coordinates.includeZ = true;
cfg.coordinates.zeroOrigin = true;

cfg.numerics = struct();
cfg.numerics.psdTol = 1e-10;

%% Step 1: first access builds the cache

fprintf('--- Step 1: first access builds cache for file 1 ---\n\n');

mode1a = loadOWNSModeCached(file1, cfg);

cacheFile1 = fullfile(cacheDir, dir(fullfile(cacheDir, '*.mat')).name);
info1a = dir(cacheFile1);

if isempty(info1a)
    error('runOWNSCacheStalenessTest:CacheNotBuilt', ...
        'No cache file was created for file 1.');
end

%% Step 2: persistence -- unchanged file must not be rebuilt

fprintf('\n--- Step 2: repeated access must reuse cache (no rebuild) ---\n\n');

mode1b = loadOWNSModeCached(file1, cfg);

info1b = dir(cacheFile1);

if info1b.datenum ~= info1a.datenum
    error('runOWNSCacheStalenessTest:UnexpectedRebuild', ...
        'Cache file for file 1 was rewritten even though the source did not change.');
end

if ~isequal(mode1a.A{1}, mode1b.A{1})
    error('runOWNSCacheStalenessTest:InconsistentReload', ...
        'Cached reload of file 1 does not match the original build.');
end

fprintf('OK: cache reused, file unchanged since Step 1.\n');

%% Step 3: staleness -- touching the source file must trigger a rebuild

fprintf('\n--- Step 3: modified source file must invalidate cache ---\n\n');

status = system(sprintf('touch -d "+1 hour" ''%s''', file1));
if status ~= 0
    error('runOWNSCacheStalenessTest:TouchFailed', ...
        'Could not update the modification time of file 1.');
end

mode1c = loadOWNSModeCached(file1, cfg);

info1c = dir(cacheFile1);

if info1c.datenum == info1b.datenum
    error('runOWNSCacheStalenessTest:StaleCacheNotRebuilt', ...
        'Cache file for file 1 was not rebuilt after the source changed.');
end

if ~isequal(mode1c.A{1}, mode1a.A{1})
    error('runOWNSCacheStalenessTest:RebuildMismatch', ...
        'Rebuilt cache for file 1 does not match the original content.');
end

fprintf('OK: stale cache detected and rebuilt after source modification.\n');

%% Step 4: a newly seen file is cached without disturbing existing entries

fprintf('\n--- Step 4: new file is cached; existing entries are untouched ---\n\n');

infoBeforeNewFile = dir(cacheFile1);

mode2 = loadOWNSModeCached(file2, cfg); %#ok<NASGU>

cacheListing = dir(fullfile(cacheDir, '*.mat'));
if numel(cacheListing) ~= 2
    error('runOWNSCacheStalenessTest:NewFileNotCached', ...
        'Expected 2 cached files after adding file 2, found %d.', ...
        numel(cacheListing));
end

infoAfterNewFile = dir(cacheFile1);

if infoAfterNewFile.datenum ~= infoBeforeNewFile.datenum
    error('runOWNSCacheStalenessTest:ExistingEntryDisturbed', ...
        'Caching a new file rebuilt an unrelated existing cache entry.');
end

fprintf('OK: file 2 cached; file 1''s cache entry was left untouched.\n');

%% Clean up

rmdir(tempDataDir, 's');
rmdir(cacheDir, 's');

fprintf('\nPASSED: runOWNSCacheStalenessTest\n');
