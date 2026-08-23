% runHermiteMaximumTests.m
%
% Verification of cubic-Hermite continuous maximum detection.
%
% MATLAB version: R2020b

clear;
%clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('Cubic-Hermite maximum regression tests\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  1. Generic configuration
% ========================================================================

cfg = struct();

cfg.gaussianConvention = 'real';
cfg.factorOrientation = 'synthesis';
cdf.synthetic.Nx = 41;

cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.maxDetection.method = 'hermite';

%% ========================================================================
%  2. Exact cubic polynomial test
% ========================================================================

fprintf('Test 1: exact cubic polynomial...\n');

% Polynomial:
%
%   a(x) = -2*x^3 + 3*x^2 + 0.2*x + 0.1
%
% on [0,1].

xLeft = 0;
xRight = 1;

a = @(x) -2*x.^3 + 3*x.^2 + 0.2*x + 0.1;
ap = @(x) -6*x.^2 + 6*x + 0.2;

valueLeft = a(xLeft);
valueRight = a(xRight);

derivativeLeft = ap(xLeft);
derivativeRight = ap(xRight);

[mHermite, xHermite] = cubicHermiteIntervalMaximum( ...
    xLeft, xRight, ...
    valueLeft, valueRight, ...
    derivativeLeft, derivativeRight);

% Exact stationary roots.
stationaryRoots = roots([-6, 6, 0.2]);
candidateX = [xLeft; xRight; ...
    stationaryRoots( ...
        isreal(stationaryRoots) ...
        & stationaryRoots > xLeft ...
        & stationaryRoots < xRight)];

candidateValues = a(candidateX);
[mExact, indexExact] = max(candidateValues);
xExact = candidateX(indexExact);

fprintf('  Hermite maximum = %.16e at x = %.16e\n', ...
    mHermite, xHermite);
fprintf('  Exact maximum   = %.16e at x = %.16e\n', ...
    mExact, xExact);

assert(abs(mHermite - mExact) < 1e-13, ...
    'Hermite interpolation did not reproduce the cubic maximum.');

assert(abs(xHermite - xExact) < 1e-13, ...
    'Hermite interpolation did not reproduce the cubic maximizer.');

fprintf('  passed.\n\n');

%% ========================================================================
%  3. Synthetic Gprime verification
% ========================================================================

fprintf('Test 2: exact synthetic Gprime...\n');

cfg.synthetic.Nx = 301;
cfg.synthetic.xMax = 20;

[x, B, H, eThresh, meta] = makeSyntheticProblem(cfg);

G = buildGFamily(B, H, eThresh, cfg);

[Gprime, derivativeInfo] = buildGprimeFamily( ...
    x, G, meta, cfg, []);

assert(derivativeInfo.usedExactDerivative, ...
    'The synthetic test should use exact derivative metadata.');

maxDerivativeError = 0;

for n = 1:numel(x)
    maxDerivativeError = max( ...
        maxDerivativeError, ...
        norm( ...
            Gprime{n} - meta.GprimeExact{n}, ...
            'fro'));
end

fprintf('  max exact derivative error = %.6e\n', ...
    maxDerivativeError);

assert(maxDerivativeError < 1e-13, ...
    'Exact synthetic Gprime data were not recovered.');

derivativeDiagnostics = ...
    validateGprimeFamily(x, G, Gprime);

fprintf('  centered-FD diagnostic     = %.6e\n', ...
    derivativeDiagnostics.maxFiniteDifferenceDefect);

fprintf('  passed.\n\n');

%% ========================================================================
%  4. Compare grid and Hermite maxima on the standard grid
% ========================================================================

fprintf('Test 3: grid versus Hermite transition CDF...\n');

angularCfg = struct();
angularCfg.method = 'circle';
angularCfg.numDirections = 4096;
angularCfg.randomSeed = 1;

r = size(B{1}, 2);

[U, weights] = generateAngularRule(r, angularCfg);

maxGrid = computeRunningMaxGrid(G, U, cfg);

maxHermite = computeRunningMaxHermite( ...
    x, G, Gprime, U, cfg);

validateRunningMax(maxGrid.m);
validateRunningMax(maxHermite.m);

probGrid = computeTransitionCDF( ...
    maxGrid.m, weights, r, cfg);

probHermite = computeTransitionCDF( ...
    maxHermite.m, weights, r, cfg);

fprintf('  grid F(xmax)    = %.12e\n', probGrid.F(end));
fprintf('  Hermite F(xmax) = %.12e\n', probHermite.F(end));
fprintf('  difference      = %.12e\n', ...
    probHermite.F(end) - probGrid.F(end));

% Since Hermite can resolve interior peaks, it should not produce a
% smaller running maximum or transition CDF than endpoint sampling.
minimumMaximumDifference = ...
    min(maxHermite.m(:) - maxGrid.m(:));

minimumCdfDifference = ...
    min(probHermite.F - probGrid.F);

fprintf('  min[mHermite - mGrid] = %.12e\n', ...
    minimumMaximumDifference);
fprintf('  min[FHermite - FGrid] = %.12e\n', ...
    minimumCdfDifference);

assert(minimumMaximumDifference > -1e-12, ...
    'Hermite running maxima are unexpectedly below grid maxima.');

assert(minimumCdfDifference > -1e-12, ...
    'Hermite transition CDF is unexpectedly below the grid CDF.');

fprintf('  passed.\n\n');

%% ========================================================================
%  5. Dense-grid reference
% ========================================================================

fprintf('Test 4: Hermite result versus dense-grid reference...\n');

cfgDense = cfg;
cfgDense.synthetic.Nx = 6001;

[xDense, BDense, HDense, eThreshDense] = ...
    makeSyntheticProblem(cfgDense);

GDense = buildGFamily( ...
    BDense, HDense, eThreshDense, cfgDense);

% Use a moderate angular resolution to keep the dense reference compact.
angularCfg.numDirections = 2048;

[UDense, weightsDense] = generateAngularRule( ...
    r, angularCfg);

maxDense = computeRunningMaxGrid( ...
    GDense, UDense, cfgDense);

probDense = computeTransitionCDF( ...
    maxDense.m, weightsDense, r, cfgDense);

% Recompute the coarse-grid Hermite result using the same angular rule.
maxHermiteDenseRule = computeRunningMaxHermite( ...
    x, G, Gprime, UDense, cfg);

probHermiteDenseRule = computeTransitionCDF( ...
    maxHermiteDenseRule.m, weightsDense, r, cfg);

terminalDifference = ...
    abs(probHermiteDenseRule.F(end) - probDense.F(end));

fprintf('  coarse Hermite F(xmax) = %.12e\n', ...
    probHermiteDenseRule.F(end));
fprintf('  dense grid F(xmax)     = %.12e\n', ...
    probDense.F(end));
fprintf('  terminal difference    = %.12e\n', ...
    terminalDifference);

% Compare the dense reference at the coarse-grid stations.
denseAtCoarse = interp1( ...
    xDense, probDense.F, x, 'linear');

completeCdfDifference = max(abs( ...
    probHermiteDenseRule.F - denseAtCoarse));

fprintf('  max complete-CDF difference = %.12e\n', ...
    completeCdfDifference);

% This tolerance can be tightened after observing the first run.
assert(terminalDifference < 1e-6, ...
    'Hermite terminal probability does not agree with the dense grid.');

assert(completeCdfDifference < 2e-5, ...
    'Hermite CDF does not agree with the dense-grid reference.');

fprintf('  passed.\n\n');

%% ========================================================================
%  6. Finite-difference Gprime fallback convergence
% ========================================================================

fprintf('Test 5: finite-difference Gprime fallback...\n');

metaWithoutDerivative = rmfield(meta, 'GprimeExact');
metaWithoutDerivative.hasGprime = false;

[GprimeFD, infoFD] = buildGprimeFamily( ...
    x, G, metaWithoutDerivative, cfg, []);

assert(infoFD.usedFiniteDifference, ...
    'Finite-difference derivative fallback was not selected.');

maxHermiteFD = computeRunningMaxHermite( ...
    x, G, GprimeFD, U, cfg);

probHermiteFD = computeTransitionCDF( ...
    maxHermiteFD.m, weights, r, cfg);

fdCdfDifference = max(abs( ...
    probHermiteFD.F - probHermite.F));

fprintf('  max exact-vs-FD Hermite CDF difference = %.12e\n', ...
    fdCdfDifference);

% The standard grid is already fairly fine. This tolerance can be
% adjusted after the first observed result.
assert(fdCdfDifference < 1e-4, ...
    'Finite-difference Gprime gives an unexpectedly different CDF.');

fprintf('  passed.\n\n');

%% ========================================================================
%  7. Diagnostic plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Hermite maximum verification');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
hold on;

plot(x, probGrid.F, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Stored stations');

plot(x, probHermite.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'Cubic Hermite');

xlabel('x');
ylabel('F_{X_{tr}}(x)');
title('Endpoint versus continuous maxima');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, probHermite.F - probGrid.F, ...
    'LineWidth', 1.5);

xlabel('x');
ylabel('\Delta F');
title('Hermite correction');
grid on;
box on;

nexttile;

directionIndex = round(size(U, 2) / 7);

plot(x, maxHermite.a(directionIndex, :), 'o-', ...
    'LineWidth', 1, ...
    'MarkerSize', 3, ...
    'DisplayName', 'Endpoint gain');

hold on;

plot( ...
    maxHermite.intervalX(directionIndex, :), ...
    maxHermite.intervalMax(directionIndex, :), ...
    '.', ...
    'MarkerSize', 10, ...
    'DisplayName', 'Interval maximum');

xlabel('x');
ylabel('a_u(x)');
title('Example directional interval maxima');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, maxGrid.m(directionIndex, :), '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Grid running maximum');

hold on;

plot(x, maxHermite.m(directionIndex, :), ...
    'LineWidth', 2, ...
    'DisplayName', 'Hermite running maximum');

xlabel('x');
ylabel('m_u(x)');
title('Example directional running maximum');
legend('Location', 'best');
grid on;
box on;

sgtitle('Continuous streamwise maximum verification');

fprintf('============================================================\n');
fprintf('All cubic-Hermite maximum tests passed.\n');
fprintf('============================================================\n');