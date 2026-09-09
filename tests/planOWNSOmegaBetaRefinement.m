% planOWNSOmegaBetaRefinement.m
%
% Build two energy-based omega-beta refinement plans from reduced caches.
% The proxy is the streamwise integral of mean modal energy, trace(A(x)).
% It is not a first-transition CDF convergence test.

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
cacheDirectory = fullfile(projectRoot, 'data', 'reducedCache');
outputDirectory = fullfile(projectRoot, 'output');

omegaLimits = [];
betaLimits = [0, 1];
numOmegaCells = 200;
numBetaCells = 200;
targetFraction = 0.99;

if ~isfolder(cacheDirectory)
    error('planOWNSOmegaBetaRefinement:CacheDirectoryNotFound', ...
        'Cache directory not found: %s', cacheDirectory);
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

cacheFiles = dir(fullfile(cacheDirectory, '*.mat'));
cacheFiles = cacheFiles(~[cacheFiles.isdir]);
[~, order] = sort({cacheFiles.name});
cacheFiles = cacheFiles(order);

M = numel(cacheFiles);
if M < 3
    error('planOWNSOmegaBetaRefinement:TooFewModes', ...
        'At least 3 cached modes are required.');
end

omega = zeros(M, 1);
beta = zeros(M, 1);
energyProxy = zeros(M, 1);
sourceFile = cell(M, 1);

for i = 1:M
    loaded = load(fullfile(cacheFiles(i).folder, cacheFiles(i).name), 'mode');
    mode = loaded.mode;
    omega(i) = mode.meta.omega;
    beta(i) = mode.meta.beta;
    sourceFile{i} = mode.meta.source;
    energyProxy(i) = trapz(mode.x(:), cellfun(@trace, mode.A));
end

if isempty(omegaLimits)
    omegaLimits = [min(omega), 0.1];
end

if any(energyProxy < 0) || any(~isfinite(energyProxy))
    error('planOWNSOmegaBetaRefinement:InvalidEnergyProxy', ...
        'Cached integrated mean energies must be finite and nonnegative.');
end

points = [omega, beta];
if size(unique(points, 'rows'), 1) ~= M
    error('planOWNSOmegaBetaRefinement:DuplicateModes', ...
        'Each cached omega-beta pair must be unique.');
end

[cellPoints, cellArea] = localCellCenters( ...
    omegaLimits, betaLimits, numOmegaCells, numBetaCells);
nearestMode = localNearestMode(cellPoints, points, omegaLimits, betaLimits);
cellEnergy = cellArea * energyProxy(nearestMode);
totalEnergy = sum(cellEnergy);

if totalEnergy <= 0
    error('planOWNSOmegaBetaRefinement:ZeroEnergyProxy', ...
        'The nearest-cell energy proxy is zero over the integration domain.');
end

% Completion plan: leave-one-out disagreement weighted by represented energy.
leaveOneOut = localLeaveOneOutRelativeError( ...
    points, energyProxy, omegaLimits, betaLimits);
cellScore = cellEnergy .* leaveOneOut(nearestMode);
completion = localCompletionPlan( ...
    cellPoints, cellEnergy, cellScore, nearestMode, omega, beta, targetFraction);

% New-data plan: greedily remove the least consequential cached mode. The
% retained modes reproduce the full-cache nearest-cell energy integral to 1%.
retained = localGreedyRetainedModes( ...
    cellPoints, cellEnergy, cellArea, energyProxy, points, ...
    omegaLimits, betaLimits, targetFraction);
newGrid = table(omega(retained), beta(retained), energyProxy(retained), ...
    sourceFile(retained), ...
    'VariableNames', {'omega', 'beta', 'integratedMeanEnergy', 'sourceFile'});

writetable(completion, fullfile(outputDirectory, ...
    'omega_beta_completion_plan.csv'));
writetable(newGrid, fullfile(outputDirectory, ...
    'omega_beta_new_grid_plan.csv'));

fprintf('Modes analyzed: %d\n', M);
fprintf('Nearest-cell proxy energy: %.6e\n', totalEnergy);
fprintf('Completion candidates through %.1f%% score: %d\n', ...
    100 * targetFraction, height(completion));
fprintf('Retained modes at %.1f%% proxy-energy agreement: %d\n', ...
    100 * targetFraction, height(newGrid));
fprintf('Wrote: %s\n', outputDirectory);

function [points, area] = localCellCenters(omegaLimits, betaLimits, nOmega, nBeta)
dOmega = diff(omegaLimits) / nOmega;
dBeta = diff(betaLimits) / nBeta;
omega = omegaLimits(1) + dOmega * ((0:(nOmega - 1)) + 0.5);
beta = betaLimits(1) + dBeta * ((0:(nBeta - 1)) + 0.5);
[omegaGrid, betaGrid] = ndgrid(omega, beta);
points = [omegaGrid(:), betaGrid(:)];
area = dOmega * dBeta;
end

function nearest = localNearestMode(query, samples, omegaLimits, betaLimits)
queryScaled = localScale(query, omegaLimits, betaLimits);
sampleScaled = localScale(samples, omegaLimits, betaLimits);
distanceSquared = (queryScaled(:, 1) - sampleScaled(:, 1).').^2 + ...
    (queryScaled(:, 2) - sampleScaled(:, 2).').^2;
[~, nearest] = min(distanceSquared, [], 2);
end

function scaled = localScale(points, omegaLimits, betaLimits)
scaled = [(points(:, 1) - omegaLimits(1)) / diff(omegaLimits), ...
    (points(:, 2) - betaLimits(1)) / diff(betaLimits)];
end

function errorEstimate = localLeaveOneOutRelativeError( ...
    points, energy, omegaLimits, betaLimits)
M = numel(energy);
scaled = localScale(points, omegaLimits, betaLimits);
distanceSquared = (scaled(:, 1) - scaled(:, 1).').^2 + ...
    (scaled(:, 2) - scaled(:, 2).').^2;
distanceSquared(1:(M + 1):end) = inf;
[~, nearest] = min(distanceSquared, [], 2);
energyFloor = max(max(energy) * 1e-12, eps);
errorEstimate = abs(energy - energy(nearest)) ./ ...
    max(energy, energyFloor);
end

function plan = localCompletionPlan( ...
    cells, cellEnergy, cellScore, nearest, omega, beta, targetFraction)
M = numel(omega);
candidateOmega = zeros(M, 1);
candidateBeta = zeros(M, 1);
representedEnergy = zeros(M, 1);
unresolvedScore = zeros(M, 1);

for i = 1:M
    owned = nearest == i;
    representedEnergy(i) = sum(cellEnergy(owned));
    unresolvedScore(i) = sum(cellScore(owned));
    [~, localIndex] = max(cellScore(owned));
    ownedIndices = find(owned);
    candidateOmega(i) = cells(ownedIndices(localIndex), 1);
    candidateBeta(i) = cells(ownedIndices(localIndex), 2);
end

[unresolvedScore, order] = sort(unresolvedScore, 'descend');
representedEnergy = representedEnergy(order);
candidateOmega = candidateOmega(order);
candidateBeta = candidateBeta(order);
omega = omega(order);
beta = beta(order);

if sum(unresolvedScore) <= 0
    plan = table([], [], [], [], [], [], [], [], ...
        'VariableNames', {'rank', 'candidateOmega', 'candidateBeta', ...
        'nearestOmega', 'nearestBeta', 'representedEnergy', ...
        'unresolvedEnergyScore', 'cumulativeScoreFraction'});
    return;
end

keep = cumsum(unresolvedScore) <= targetFraction * sum(unresolvedScore);
if any(~keep) && ~isempty(keep)
    keep(find(~keep, 1, 'first')) = true;
end

plan = table((1:nnz(keep)).', candidateOmega(keep), candidateBeta(keep), ...
    omega(keep), beta(keep), representedEnergy(keep), unresolvedScore(keep), ...
    cumsum(unresolvedScore(keep)) / max(sum(unresolvedScore), eps), ...
    'VariableNames', {'rank', 'candidateOmega', 'candidateBeta', ...
    'nearestOmega', 'nearestBeta', 'representedEnergy', ...
    'unresolvedEnergyScore', 'cumulativeScoreFraction'});
end

function retained = localGreedyRetainedModes( ...
    cells, cellEnergy, cellArea, energy, samples, omegaLimits, ...
    betaLimits, targetFraction)
M = size(samples, 1);
retained = true(M, 1);
referenceEnergy = sum(cellEnergy);

while nnz(retained) > 1
    active = find(retained);
    bestIndex = 0;
    bestDefect = inf;

    for j = 1:numel(active)
        trial = retained;
        trial(active(j)) = false;
        nearest = localNearestMode(cells, samples(trial, :), ...
            omegaLimits, betaLimits);
        trialEnergy = sum(cellArea * energy(trial(nearest)));
        defect = abs(trialEnergy - referenceEnergy) / referenceEnergy;
        if defect < bestDefect
            bestDefect = defect;
            bestIndex = active(j);
        end
    end

    if bestDefect > 1 - targetFraction
        break;
    end
    retained(bestIndex) = false;
end

end