% runSyntheticRegressionTests.m
%
% Regression and convergence tests for the rotating two-dimensional
% synthetic first-transition problem.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'examples')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('Synthetic first-transition regression tests\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  1. Configuration
% ========================================================================

cfg = struct();

cfg.gaussianConvention = 'real';
cfg.factorOrientation = 'synthesis';
cfg.isPhysicalStateComplex = true;

cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.maxDetection.method = 'grid';

%% ========================================================================
%  2. Generate the test problem
% ========================================================================

[x, B, H, eThresh, meta] = makeSyntheticProblem(cfg);

B = convertFactorOrientation(B, cfg.factorOrientation);

validateProblemData(x, B, H, eThresh, cfg);

G = buildGFamily(B, H, eThresh, cfg);

gDiagnostics = validateGFamily(G, cfg);

Nx = numel(x);
r = size(B{1}, 2);

assert(r == 2, ...
    'The current synthetic regression problem is expected to have r=2.');

fprintf('Generated synthetic problem:\n');
fprintf('  Nx = %d\n', Nx);
fprintf('  r  = %d\n\n', r);

%% ========================================================================
%  3. Check reconstruction of the analytical synthetic G(x)
% ========================================================================

gError = zeros(Nx, 1);

for n = 1:Nx
    gError(n) = norm( ...
        G{n} - meta.GExact{n}, ...
        'fro');
end

maxGError = max(gError);

fprintf('G reconstruction check:\n');
fprintf('  max ||G - GExact||_F = %.6e\n\n', maxGError);

assert(maxGError < 1e-12, ...
    'The constructed G matrices do not match the synthetic reference.');

%% ========================================================================
%  4. Circular angular-convergence study
% ========================================================================

KValues = [64, 128, 256, 512, 1024, 2048, 4096, 8192];

numK = numel(KValues);

terminalProbability = zeros(numK, 1);
cdfDifference = nan(numK, 1);

cdfPrevious = [];

fprintf('Circular angular-convergence study:\n');
fprintf('%10s %20s %20s\n', ...
    'K', 'F(xmax)', 'max CDF change');

for j = 1:numK
    angularCfg = struct();
    angularCfg.method = 'circle';
    angularCfg.numDirections = KValues(j);
    angularCfg.randomSeed = 1;

    [U, weights] = generateAngularRule(r, angularCfg);
    validateAngularRule(U, weights);

    maxData = computeRunningMaxGrid(G, U, cfg);
    validateRunningMax(maxData.m);

    prob = computeTransitionCDF( ...
        maxData.m, weights, r, cfg);

    validateProbabilityCurves(prob, cfg);

    terminalProbability(j) = prob.F(end);

    if ~isempty(cdfPrevious)
        cdfDifference(j) = max(abs(prob.F - cdfPrevious));
    end

    fprintf('%10d %20.12e %20.12e\n', ...
        KValues(j), ...
        terminalProbability(j), ...
        cdfDifference(j));

    cdfPrevious = prob.F;
end

fprintf('\n');

FReference = terminalProbability(end);
angularDifference = cdfDifference(end);

fprintf('Highest-resolution circular result:\n');
fprintf('  K                 = %d\n', KValues(end));
fprintf('  F(xmax)           = %.12e\n', FReference);
fprintf('  max CDF refinement change = %.12e\n\n', ...
    angularDifference);

% The circular trapezoidal rule should converge rapidly for this smooth
% low-dimensional test. This tolerance can be tightened after observing
% the first convergence table.
assert(angularDifference < 1e-8, ...
    'Circular angular quadrature has not converged to the test tolerance.');

%% ========================================================================
%  5. Independent trajectory Monte Carlo check
% ========================================================================

numSamples = 2e5;
mcSeed = 24681357;
mcBlockSize = 20000;

fprintf('Running direct trajectory Monte Carlo:\n');
fprintf('  samples = %d\n', numSamples);

mc = computeTrajectoryMonteCarlo( ...
    G, numSamples, mcSeed, mcBlockSize);

fprintf('  Monte Carlo F(xmax) = %.12e\n', mc.F(end));
fprintf('  Circular     F(xmax) = %.12e\n', FReference);
fprintf('  Difference           = %.12e\n', ...
    abs(mc.F(end) - FReference));
fprintf('  MC standard error    = %.12e\n\n', ...
    mc.standardError(end));

% Use the highest-resolution circular CDF as the deterministic reference.
angularCfg.method = 'circle';
angularCfg.numDirections = KValues(end);

[UReference, weightsReference] = ...
    generateAngularRule(r, angularCfg);

maxReference = computeRunningMaxGrid( ...
    G, UReference, cfg);

probReference = computeTransitionCDF( ...
    maxReference.m, weightsReference, r, cfg);

difference = mc.F - probReference.F;

% Pointwise standardized discrepancy. Protect zero standard errors.
standardErrorFloor = 1 / numSamples;
standardizedDifference = ...
    abs(difference) ./ max(mc.standardError, standardErrorFloor);

maxStandardizedDifference = max(standardizedDifference);

fprintf('Monte Carlo comparison over complete CDF:\n');
fprintf('  max absolute CDF difference = %.12e\n', ...
    max(abs(difference)));
fprintf('  max standardized difference = %.6f standard errors\n\n', ...
    maxStandardizedDifference);

% A pointwise 5-sigma bound is intentionally loose because the maximum is
% taken over many correlated streamwise stations.
assert(maxStandardizedDifference < 5, ...
    'Monte Carlo and radial-directional CDFs disagree unexpectedly.');

%% ========================================================================
%  6. Local exceedance consistency
% ========================================================================

pLocal = computeLocalExceedance( ...
    G, UReference, weightsReference, r, cfg);

localViolation = max(pLocal - probReference.F);

fprintf('Local exceedance consistency:\n');
fprintf('  max[pLocal - F_Xtr] = %.12e\n\n', localViolation);

assert(localViolation < 1e-12, ...
    'Local exceedance probability exceeds first-transition CDF.');

assert(abs(pLocal(1) - probReference.F(1)) < 1e-12, ...
    'Local and first-transition probabilities must agree at the inlet.');

%% ========================================================================
%  7. Probability conservation
% ========================================================================

totalProbability = ...
    sum(probReference.pInterval) ...
    + probReference.pCensored;

fprintf('Probability conservation:\n');
fprintf('  total probability = %.16f\n', totalProbability);
fprintf('  defect            = %.6e\n\n', ...
    abs(totalProbability - 1));

assert(abs(totalProbability - 1) < 1e-12, ...
    'Probability is not conserved.');

%% ========================================================================
%  8. Regression target
% ========================================================================

% This value is based on the current synthetic problem, the current
% 301-point streamwise grid, endpoint-only maxima, and a highly resolved
% circular angular rule.
%
% If the synthetic problem or continuous maximum treatment is changed,
% update this reference only after independently verifying the new result.

regressionTarget = 3.461435511754e-01;
regressionTolerance = 5e-7;

regressionError = abs( ...
    probReference.F(end) - regressionTarget);

fprintf('Terminal-probability regression check:\n');
fprintf('  computed = %.12e\n', probReference.F(end));
fprintf('  target   = %.12e\n', regressionTarget);
fprintf('  error    = %.12e\n\n', regressionError);

assert(regressionError < regressionTolerance, ...
    'Terminal transition probability changed from its regression target.');

%% ========================================================================
%  9. Diagnostic plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'Synthetic regression diagnostics');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

semilogx(KValues, terminalProbability, 'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Number of angular points K');
ylabel('F_{X_{tr}}(x_{max})');
title('Angular convergence');
grid on;
box on;

nexttile;

semilogy(KValues(2:end), cdfDifference(2:end), 'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Number of angular points K');
ylabel('Maximum CDF change');
title('Angular refinement difference');
grid on;
box on;

nexttile;
hold on;

plot(x, probReference.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'Radial--directional');

plot(x, mc.F, '--', ...
    'LineWidth', 1.25, ...
    'DisplayName', 'Trajectory Monte Carlo');

xlabel('x');
ylabel('F_{X_{tr}}(x)');
title('Independent probability comparison');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, pLocal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Local exceedance');

hold on;

plot(x, probReference.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'First-transition CDF');

xlabel('x');
ylabel('Probability');
title('Effect of streamwise memory');
legend('Location', 'best');
grid on;
box on;

sgtitle('Synthetic first-transition regression tests');

fprintf('============================================================\n');
fprintf('All synthetic regression tests passed.\n');
fprintf('============================================================\n');