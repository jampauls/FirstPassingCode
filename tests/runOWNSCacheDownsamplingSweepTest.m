% runOWNSCacheDownsamplingSweepTest.m
%
% Sweeps the tunable streamwise-station downsampling parameters used when
% building the reduced cache (cfg.owns.numCacheStations and
% cfg.owns.cacheXCutoffCoefficient, see downsampleOWNSStations) across 5
% resolutions:
%
%   - the low-resolution (few-station) sweep points span the entire
%     streamwise grid;
%   - the high-resolution (many-station) sweep points instead restrict
%     the eligible domain to a 1/sqrt(omega) cutoff, so the extra station
%     budget is concentrated near the inlet rather than spread across a
%     domain most high-frequency modes never reach.
%
% Uses two of the existing OWNS data files with different omega so the
% omega-dependent cutoff can be checked directly: the higher-omega file
% must retain a smaller streamwise span than the lower-omega file once
% the cutoff is enabled.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('OWNS cache downsampling sweep test\n');
fprintf('============================================================\n\n');

dataDir = fullfile(projectRoot, 'data');
cacheDir = fullfile(tempdir, 'FirstPassingCodeDownsamplingSweepCache');

if isfolder(cacheDir)
    rmdir(cacheDir, 's');
end

baseCfg = struct();
baseCfg.owns = struct();
baseCfg.owns.reducedCacheDirectory = cacheDir;
baseCfg.owns.solutionVariable = 'solution';
baseCfg.owns.numEnergyVariables = 5;
baseCfg.owns.numStateVariables = 6;
baseCfg.owns.validateWeights = true;
baseCfg.owns.weightTolerance = 1e-12;
baseCfg.owns.verbose = false;

baseCfg.coordinates = struct();
baseCfg.coordinates.method = 'wallArcLength';
baseCfg.coordinates.referenceIndex = 1;
baseCfg.coordinates.includeZ = true;
baseCfg.coordinates.zeroOrigin = true;

baseCfg.numerics = struct();
baseCfg.numerics.psdTol = 1e-10;

fileList = discoverOWNSEnsembleFiles(dataDir);
if numel(fileList) < 2
    error('runOWNSCacheDownsamplingSweepTest:TooFewFiles', ...
        'Expected at least 2 OWNS mode files in data/, found %d.', ...
        numel(fileList));
end

%% Probe every file cheaply (2 stations, no cutoff) for omega and x-span

fprintf('Probing available mode files for omega and streamwise span...\n\n');

probeCfg = baseCfg;
probeCfg.owns.numCacheStations = 2;
probeCfg.owns.cacheXCutoffCoefficient = [];

omegaByFile = zeros(numel(fileList), 1);
fullSpanByFile = zeros(numel(fileList), 1);

for k = 1:numel(fileList)
    probeMode = loadOWNSModeCached(fileList{k}, probeCfg);
    omegaByFile(k) = probeMode.meta.omega;
    fullSpanByFile(k) = probeMode.x(end) - probeMode.x(1);
    [~, fname] = fileparts(fileList{k});
    fprintf('  %-28s omega = %10.4g   full span = %.6g\n', ...
        fname, omegaByFile(k), fullSpanByFile(k));
end
fprintf('\n');

rmdir(cacheDir, 's'); % discard probe cache; the sweep below rebuilds it

[~, loFile] = min(abs(omegaByFile));
[~, hiFile] = max(abs(omegaByFile));

if loFile == hiFile
    error('runOWNSCacheDownsamplingSweepTest:NoOmegaContrast', ...
        'All discovered mode files share the same omega; cannot test the cutoff.');
end

fileLow = fileList{loFile};
fileHigh = fileList{hiFile};

fprintf('Low-omega contrast file:  %s (omega = %.4g)\n', fileLow, omegaByFile(loFile));
fprintf('High-omega contrast file: %s (omega = %.4g)\n\n', fileHigh, omegaByFile(hiFile));

% Coefficient chosen relative to the observed problem scales so the
% cutoff meaningfully truncates the domain for a typical mode, without
% depending on any hard-coded absolute length scale.
typicalScale = median(fullSpanByFile .* sqrt(abs(omegaByFile)));
cutoffCoefficient = 0.4 * typicalScale;

%% Sweep configuration: 5 station counts, cutoff enabled only for the top 2

sweepNumStations = [20, 40, 80, 160, 320];
sweepUseCutoff   = [false, false, false, true, true];

fprintf('Sweep station counts: %s\n', mat2str(sweepNumStations));
fprintf('Cutoff enabled for  : %s\n\n', mat2str(sweepUseCutoff));

for i = 1:numel(sweepNumStations)
    numStations = sweepNumStations(i);
    useCutoff = sweepUseCutoff(i);

    cfg = baseCfg;
    cfg.owns.numCacheStations = numStations;

    if useCutoff
        cfg.owns.cacheXCutoffCoefficient = cutoffCoefficient;
    else
        cfg.owns.cacheXCutoffCoefficient = [];
    end

    fprintf('--- Sweep point %d: numCacheStations = %d, cutoff = %d ---\n', ...
        i, numStations, useCutoff);

    modeLow = loadOWNSModeCached(fileLow, cfg);
    modeHigh = loadOWNSModeCached(fileHigh, cfg);

    %% Retained station count must match the requested budget (or fewer,
    %% only if the eligible/raw grid itself has fewer stations)

    if numel(modeLow.x) > numStations || numel(modeHigh.x) > numStations
        error('runOWNSCacheDownsamplingSweepTest:TooManyStations', ...
            'Sweep point %d retained more stations than requested.', i);
    end

    if numel(modeLow.x) < 2 || numel(modeHigh.x) < 2
        error('runOWNSCacheDownsamplingSweepTest:TooFewStations', ...
            'Sweep point %d retained fewer than 2 stations.', i);
    end

    %% Streamwise coordinate must be strictly increasing (already
    %% enforced inside extractOWNSCoordinate; re-check defensively here)

    if any(diff(modeLow.x) <= 0) || any(diff(modeHigh.x) <= 0)
        error('runOWNSCacheDownsamplingSweepTest:NonIncreasing', ...
            'Sweep point %d produced a non-increasing streamwise grid.', i);
    end

    spanLow = modeLow.x(end) - modeLow.x(1);
    spanHigh = modeHigh.x(end) - modeHigh.x(1);

    if ~useCutoff
        %% Low-resolution sweep points must span the entire raw grid
        fullSpanTol = 1e-9 * max(fullSpanByFile(loFile), 1);
        if abs(spanLow - fullSpanByFile(loFile)) > fullSpanTol
            error('runOWNSCacheDownsamplingSweepTest:NotFullSpan', ...
                'Sweep point %d (no cutoff) did not span the full grid for %s.', ...
                i, fileLow);
        end

        fullSpanTol = 1e-9 * max(fullSpanByFile(hiFile), 1);
        if abs(spanHigh - fullSpanByFile(hiFile)) > fullSpanTol
            error('runOWNSCacheDownsamplingSweepTest:NotFullSpan', ...
                'Sweep point %d (no cutoff) did not span the full grid for %s.', ...
                i, fileHigh);
        end

        fprintf('  OK: both modes span the full streamwise grid.\n');
    else
        %% Cutoff sweep points must be strictly truncated relative to the
        %% full domain, and the higher-omega mode must be truncated more
        if spanLow >= fullSpanByFile(loFile) - eps(fullSpanByFile(loFile))
            error('runOWNSCacheDownsamplingSweepTest:CutoffNotApplied', ...
                'Sweep point %d did not truncate the domain for %s.', ...
                i, fileLow);
        end

        if spanHigh >= fullSpanByFile(hiFile) - eps(fullSpanByFile(hiFile))
            error('runOWNSCacheDownsamplingSweepTest:CutoffNotApplied', ...
                'Sweep point %d did not truncate the domain for %s.', ...
                i, fileHigh);
        end

        if spanHigh >= spanLow
            error('runOWNSCacheDownsamplingSweepTest:CutoffNotOmegaScaled', ...
                ['Sweep point %d: the higher-omega mode retained a span ', ...
                 '(%.6g) at least as large as the lower-omega mode (%.6g).'], ...
                i, spanHigh, spanLow);
        end

        fprintf(['  OK: domain truncated for both modes; higher-omega span ', ...
            '(%.6g) < lower-omega span (%.6g).\n'], spanHigh, spanLow);
    end

    fprintf('  stations retained: %s = %d, %s = %d\n\n', ...
        fileLow, numel(modeLow.x), fileHigh, numel(modeHigh.x));
end

%% Every sweep point must have produced its own distinct cache entry

cacheListing = dir(fullfile(cacheDir, '*.mat'));
expectedCacheFiles = 2 * numel(sweepNumStations);

if numel(cacheListing) ~= expectedCacheFiles
    error('runOWNSCacheDownsamplingSweepTest:CacheCollision', ...
        'Expected %d distinct cache files across the sweep, found %d.', ...
        expectedCacheFiles, numel(cacheListing));
end

fprintf('OK: %d distinct cache entries for %d sweep points x 2 files.\n\n', ...
    numel(cacheListing), numel(sweepNumStations));

%% Clean up

rmdir(cacheDir, 's');

fprintf('PASSED: runOWNSCacheDownsamplingSweepTest\n');
