% runProperComplexRegressionTests.m
%
% Analytical and trajectory-Monte-Carlo verification of the
% proper-complex first-transition implementation.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('Proper-complex regression tests\n');
fprintf('============================================================\n\n');

%% ========================================================================
% 1. Scalar complex analytical test
% ========================================================================

fprintf('Test 1: scalar proper-complex analytical case...\n');

x = linspace(0, 10, 101).';

g = ...
    0.1 ...
    + 0.9 * exp(-0.5 * ((x - 4.5) / 1.1).^2) ...
    + 0.4 * exp(-0.5 * ((x - 8.0) / 0.8).^2);

G = cell(numel(x), 1);

for n = 1:numel(x)
    G{n} = g(n);
end

% In complex dimension one, the complex unit direction is only a phase.
% The directional quadratic form is independent of that phase.
mExact = cummax(g);

% Since |z|^2 ~ Exponential(1):
%
%   F(x) = exp[-1/m(x)].
FExact = exp(-1 ./ mExact);

angularCfg = struct();
angularCfg.numDirections = 256;
angularCfg.numReplicates = 4;
angularCfg.randomSeed = 123;
angularCfg.sobolSkip = 128;
angularCfg.sobolLeap = 0;
angularCfg.scramble = true;
angularCfg.probabilityClip = 1e-12;

maxCfg = struct();
maxCfg.method = 'grid';

numericsCfg = struct();
numericsCfg.clipSmallNegativeGains = true;
numericsCfg.psdTol = 1e-12;

result = computeProperComplexRQMCCDF( ...
    x, G, [], 1, ...
    angularCfg, maxCfg, numericsCfg);

scalarError = max(abs(result.F - FExact));

fprintf('  maximum analytical CDF error = %.12e\n', ...
    scalarError);

assert(scalarError < 1e-12, ...
    'Scalar proper-complex analytical test failed.');

fprintf('  passed.\n\n');

%% ========================================================================
% 2. Fixed Hermitian matrix one-point moment test
% ========================================================================

fprintf('Test 2: proper-complex quadratic-form moments...\n');

rng(2468, 'twister');

r = 4;

Qraw = randn(r) + 1i * randn(r);
[Q, ~] = qr(Qraw);

lambda = [1.2, 0.7, 0.25, 0.05].';

A = Q * diag(lambda) * Q';

meanExact = trace(A);
varianceExact = real(trace(A * A));

numSamples = 5e5;
blockSize = 25000;
numBlocks = ceil(numSamples / blockSize);

sumEnergy = 0;
sumEnergySquared = 0;
numGenerated = 0;

for b = 1:numBlocks
    currentSize = min( ...
        blockSize, numSamples - numGenerated);

    Z = ( ...
        randn(r, currentSize) ...
        + 1i * randn(r, currentSize)) ...
        / sqrt(2);

    energy = real(sum(conj(Z) .* (A * Z), 1));

    sumEnergy = sumEnergy + sum(energy);
    sumEnergySquared = ...
        sumEnergySquared + sum(energy.^2);

    numGenerated = numGenerated + currentSize;
end

meanSample = sumEnergy / numSamples;

varianceSample = ...
    sumEnergySquared / numSamples ...
    - meanSample^2;

meanRelativeError = ...
    abs(meanSample - meanExact) / meanExact;

varianceRelativeError = ...
    abs(varianceSample - varianceExact) / varianceExact;

fprintf('  mean exact/sample     = %.8e / %.8e\n', ...
    meanExact, meanSample);
fprintf('  variance exact/sample = %.8e / %.8e\n', ...
    varianceExact, varianceSample);
fprintf('  mean relative error   = %.6e\n', ...
    meanRelativeError);
fprintf('  variance relative error = %.6e\n', ...
    varianceRelativeError);

assert(meanRelativeError < 5e-3, ...
    'Proper-complex sample mean is inconsistent.');

assert(varianceRelativeError < 1e-2, ...
    'Proper-complex sample variance is inconsistent.');

fprintf('  passed.\n\n');

%% ========================================================================
% 3. Rotating two-dimensional Hermitian process
% ========================================================================

fprintf('Test 3: RQMC versus direct proper-complex trajectories...\n');

Nx = 151;
x = linspace(0, 12, Nx).';

G = cell(Nx, 1);

for n = 1:Nx
    xn = x(n);

    lambda1 = ...
        0.05 ...
        + 0.9 * exp(-0.5 * ((xn - 4.0) / 1.0)^2);

    lambda2 = ...
        0.03 ...
        + 0.6 * exp(-0.5 * ((xn - 8.0) / 1.2)^2);

    theta = 0.18 * xn;

    % Include a genuinely complex relative phase.
    phase = 0.12 * xn;

    Rreal = [ ...
        cos(theta), -sin(theta); ...
        sin(theta),  cos(theta)];

    P = diag([1, exp(1i * phase)]);

    Ux = P * Rreal;

    G{n} = ...
        Ux * diag([lambda1, lambda2]) * Ux';

    G{n} = 0.5 * (G{n} + G{n}');
end

angularCfg.numDirections = 8192;
angularCfg.numReplicates = 8;
angularCfg.randomSeed = 9876;

rqmc = computeProperComplexRQMCCDF( ...
    x, G, [], 2, ...
    angularCfg, maxCfg, numericsCfg);

mc = computeProperComplexTrajectoryMonteCarlo( ...
    G, 3e5, 13579, 25000);

terminalDifference = ...
    abs(rqmc.F(end) - mc.F(end));

combinedScale = sqrt( ...
    rqmc.terminalStandardError^2 ...
    + mc.standardError(end)^2);

fprintf('  RQMC F(xmax) = %.12e\n', rqmc.F(end));
fprintf('  MC   F(xmax) = %.12e\n', mc.F(end));
fprintf('  difference   = %.12e\n', terminalDifference);
fprintf('  combined uncertainty scale = %.12e\n', ...
    combinedScale);

assert(terminalDifference < 5 * combinedScale, ...
    'Proper-complex RQMC disagrees with trajectory Monte Carlo.');

completeDifference = max(abs(rqmc.F - mc.F));
completeScale = max( ...
    sqrt(rqmc.standardError.^2 + mc.standardError.^2));

fprintf('  maximum complete-CDF difference = %.12e\n', ...
    completeDifference);

assert(completeDifference < 5 * completeScale, ...
    'Proper-complex complete CDF disagrees with Monte Carlo.');

fprintf('  passed.\n\n');

fprintf('============================================================\n');
fprintf('All proper-complex regression tests passed.\n');
fprintf('============================================================\n');
