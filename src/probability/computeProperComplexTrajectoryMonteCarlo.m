function mc = computeProperComplexTrajectoryMonteCarlo( ...
    G, numSamples, seed, blockSize)
%COMPUTEPROPERCOMPLEXTRAJECTORYMONTECARLO Direct complex Monte Carlo.
%
% Each realization is:
%
%   z = (x + 1i*y)/sqrt(2),
%
% where x,y are independent N(0,I).
%
% The normalized energy is:
%
%   e_n = z' * G{n} * z.
%
% MATLAB version: R2020b

if nargin < 3 || isempty(seed)
    seed = 1;
end

if nargin < 4 || isempty(blockSize)
    blockSize = 10000;
end

validateattributes(numSamples, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'numSamples');

validateattributes(blockSize, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'blockSize');

Nx = numel(G);

if Nx == 0
    error('computeProperComplexTrajectoryMonteCarlo:EmptyG', ...
        'G must contain at least one matrix.');
end

r = size(G{1}, 1);

rng(seed, 'twister');

firstTransitionIndex = ...
    repmat(Nx + 1, numSamples, 1);

numBlocks = ceil(numSamples / blockSize);

sampleMeanByStation = zeros(Nx, 1);
sampleSecondMomentByStation = zeros(Nx, 1);

for b = 1:numBlocks
    firstSample = (b - 1) * blockSize + 1;
    lastSample = min(b * blockSize, numSamples);
    blockIndices = firstSample:lastSample;
    currentBlockSize = numel(blockIndices);

    Z = ( ...
        randn(r, currentBlockSize) ...
        + 1i * randn(r, currentBlockSize)) ...
        / sqrt(2);

    transitioned = false(1, currentBlockSize);

    for n = 1:Nx
        GZ = G{n} * Z;

        energy = real(sum(conj(Z) .* GZ, 1));

        sampleMeanByStation(n) = ...
            sampleMeanByStation(n) + sum(energy);

        sampleSecondMomentByStation(n) = ...
            sampleSecondMomentByStation(n) ...
            + sum(energy.^2);

        newlyTransitioned = ...
            ~transitioned & (energy >= 1);

        if any(newlyTransitioned)
            globalIndices = blockIndices(newlyTransitioned);
            firstTransitionIndex(globalIndices) = n;
        end

        transitioned = transitioned | newlyTransitioned;

        % Do not stop early because all-station moments are also collected.
    end
end

F = zeros(Nx, 1);

for n = 1:Nx
    F(n) = mean(firstTransitionIndex <= n);
end

S = 1 - F;

sampleMeanByStation = ...
    sampleMeanByStation / numSamples;

sampleSecondMomentByStation = ...
    sampleSecondMomentByStation / numSamples;

sampleVarianceByStation = ...
    sampleSecondMomentByStation ...
    - sampleMeanByStation.^2;

standardError = sqrt( ...
    max(F .* (1 - F), 0) / numSamples);

mc = struct();
mc.F = F;
mc.S = S;
mc.standardError = standardError;
mc.firstTransitionIndex = firstTransitionIndex;
mc.pCensored = mean(firstTransitionIndex == Nx + 1);

mc.sampleMeanEnergy = sampleMeanByStation;
mc.sampleVarianceEnergy = sampleVarianceByStation;

mc.numSamples = numSamples;
mc.seed = seed;
mc.blockSize = blockSize;
mc.complexDimension = r;

end
