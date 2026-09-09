% inspectOWNSCachedSpectralGrid.m
%
% Write omega-beta coordinates from an existing reduced-cache directory.
% This script reads cached mode metadata only. It does not access raw OWNS
% solution files or create cache entries.

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

% Set this to the reduced-cache directory on the remote machine.
cacheDirectory = fullfile(projectRoot, 'data', 'reducedCache');
outputFile = fullfile(projectRoot, 'omega_beta_grid.csv');

if ~isfolder(cacheDirectory)
    error('inspectOWNSCachedSpectralGrid:CacheDirectoryNotFound', ...
        'Cache directory not found: %s', cacheDirectory);
end

cacheFiles = dir(fullfile(cacheDirectory, '*.mat'));
cacheFiles = cacheFiles(~[cacheFiles.isdir]);

if isempty(cacheFiles)
    error('inspectOWNSCachedSpectralGrid:NoCacheFiles', ...
        'No MAT-files found in cache directory: %s', cacheDirectory);
end

[~, order] = sort({cacheFiles.name});
cacheFiles = cacheFiles(order);

numFiles = numel(cacheFiles);
sourceFile = cell(numFiles, 1);
omega = zeros(numFiles, 1);
beta = zeros(numFiles, 1);

for k = 1:numFiles
    loaded = load(fullfile(cacheFiles(k).folder, cacheFiles(k).name), 'mode');

    if ~isfield(loaded, 'mode') || ~isfield(loaded.mode, 'meta') || ...
            ~isfield(loaded.mode.meta, 'omega') || ...
            ~isfield(loaded.mode.meta, 'beta')
        error('inspectOWNSCachedSpectralGrid:InvalidCacheFile', ...
            'Cache file lacks mode.meta.omega/beta: %s', cacheFiles(k).name);
    end

    sourceFile{k} = loaded.mode.meta.source;
    omega(k) = loaded.mode.meta.omega;
    beta(k) = loaded.mode.meta.beta;
end

gridTable = table(sourceFile, omega, beta, ...
    'VariableNames', {'sourceFile', 'omega', 'beta'});
writetable(gridTable, outputFile);

uniqueOmega = unique(omega);
uniqueBeta = unique(beta);
pairs = [omega, beta];
numUniquePairs = size(unique(pairs, 'rows'), 1);
expectedPairs = numel(uniqueOmega) * numel(uniqueBeta);

fprintf('Cached modes: %d\n', numFiles);
fprintf('Unique omega values: %d\n', numel(uniqueOmega));
fprintf('Unique beta values: %d\n', numel(uniqueBeta));
fprintf('Unique omega-beta pairs: %d of %d expected\n', ...
    numUniquePairs, expectedPairs);
fprintf('Complete tensor-product grid: %s\n', ...
    mat2str(numFiles == expectedPairs && numUniquePairs == expectedPairs));
fprintf('Wrote: %s\n', outputFile);