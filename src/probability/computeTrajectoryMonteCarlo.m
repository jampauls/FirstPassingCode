function mc = computeTrajectoryMonteCarlo(G, numSamples, seed, blockSize)
%COMPUTETRAJECTORYMONTECARLO Direct Monte Carlo first-transition reference.
%
%   mc = computeTrajectoryMonteCarlo(G, numSamples, seed, blockSize)
%
% For each realization:
%
%   w ~ N(0,I_r)
%   y_n = w' G(x_n) w
%
% The first-transition event by station n is:
%
%   max_{j<=n} y_j >= 1.
%
% Inputs:
%   G           cell array of r-by-r normalized energy matrices
%   numSamples  number of independent Gaussian realizations
%   seed        random-number seed
%   blockSize   samples processed per block
%
% Output:
%   mc.F                    empirical first-transition CDF
%   mc.S                    empirical survival probability
%   mc.standardError        pointwise binomial standard error
%   mc.firstTransitionIndex first-transition index for each realization
%   mc.pCensored            empirical censoring probability
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
    error('computeTrajectoryMonteCarlo:EmptyG', ...
        'G must contain at least one matrix.');
end

r = size(G{1}, 1);

for n = 1:Nx
    if any(size(G{n}) ~= [r, r])
        error('computeTrajectoryMonteCarlo:DimensionMismatch', ...
            'All G matrices must have size r-by-r.');
    end
end

rng(seed, 'twister');

% A value of Nx + 1 denotes right censoring.
firstTransitionIndex = repmat(Nx + 1, numSamples, 1);

numBlocks = ceil(numSamples / blockSize);

for b = 1:numBlocks
    firstSample = (b - 1) * blockSize + 1;
    lastSample = min(b * blockSize, numSamples);
    blockIndices = firstSample:lastSample;
    currentBlockSize = numel(blockIndices);

    % Each column is one Gaussian realization.
    W = randn(r, currentBlockSize);

    transitioned = false(1, currentBlockSize);

    for n = 1:Nx
        GW = G{n} * W;

        % Columnwise quadratic form w'G(x_n)w.
        energy = real(sum(W .* GW, 1));

        newlyTransitioned = ~transitioned & (energy >= 1);

        if any(newlyTransitioned)
            globalIndices = blockIndices(newlyTransitioned);
            firstTransitionIndex(globalIndices) = n;
        end

        transitioned = transitioned | newlyTransitioned;

        if all(transitioned)
            break;
        end
    end
end

F = zeros(Nx, 1);

for n = 1:Nx
    F(n) = mean(firstTransitionIndex <= n);
end

S = 1 - F;

standardError = sqrt( ...
    max(F .* (1 - F), 0) / numSamples);

mc = struct();
mc.F = F;
mc.S = S;
mc.standardError = standardError;
mc.firstTransitionIndex = firstTransitionIndex;
mc.pCensored = mean(firstTransitionIndex == Nx + 1);

mc.numSamples = numSamples;
mc.seed = seed;
mc.blockSize = blockSize;
mc.radialDimension = r;

end