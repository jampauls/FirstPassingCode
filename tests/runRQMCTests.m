% runRQMCTests.m
%
% Verification of scrambled-Sobol' angular integration.
%
% MATLAB version: R2020b

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'tests')));

fprintf('============================================================\n');
fprintf('RQMC angular integration tests\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  1. Problem configuration
% ========================================================================

cfg = struct();

cfg.gaussianConvention = 'real';
cfg.factorOrientation = 'synthesis';

cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

cfg.maxDetection.method = 'hermite';
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;

cfg.synthetic.Nx = 301;
cfg.synthetic.xMax = 20;

[x, B, H, eThresh, meta] = makeSyntheticProblem(cfg);

G = buildGFamily(B, H, eThresh, cfg);

[Gprime, ~] = buildGprimeFamily( ...
    x, G, meta, cfg, []);

r = size(B{1}, 2);

%% ========================================================================
%  2. Deterministic circular reference
% ========================================================================

circleCfg = struct();
circleCfg.method = 'circle';
circleCfg.numDirections = 16384;
circleCfg.randomSeed = 1;

[UCircle, wCircle] = ...
    generateAngularRule(r, circleCfg);

maxCircle = computeRunningMaxHermite( ...
    x, G, Gprime, UCircle, cfg);

probCircle = computeTransitionCDF( ...
    maxCircle.m, wCircle, r, cfg);

fprintf('Deterministic circular reference:\n');
fprintf('  F(xmax) = %.12e\n\n', probCircle.F(end));

%% ========================================================================
%  3. RQMC resolution study
% ========================================================================

KValues = [128, 256, 512, 1024, 2048, 4096];
J = 8;

terminalMean = zeros(numel(KValues), 1);
terminalSE = zeros(numel(KValues), 1);
maxCdfError = zeros(numel(KValues), 1);

fprintf('RQMC convergence study:\n');
fprintf('%8s %18s %18s %18s\n', ...
    'K', 'F(xmax)', 'terminal SE', 'max CDF error');

for j = 1:numel(KValues)
    angularCfg = struct();

    angularCfg.method = 'rqmc';
    angularCfg.numDirections = KValues(j);
    angularCfg.numReplicates = J;

    angularCfg.randomSeed = 12345;
    angularCfg.sobolSkip = 1024;
    angularCfg.sobolLeap = 0;
    angularCfg.scramble = true;
    angularCfg.probabilityClip = 1e-12;

    rqmc = computeRQMCCDF( ...
        x, G, Gprime, r, ...
        angularCfg, ...
        cfg.maxDetection, ...
        cfg.numerics);

    terminalMean(j) = rqmc.F(end);
    terminalSE(j) = rqmc.terminalStandardError;
    maxCdfError(j) = max(abs(rqmc.F - probCircle.F));

    fprintf('%8d %18.10e %18.10e %18.10e\n', ...
        KValues(j), ...
        terminalMean(j), ...
        terminalSE(j), ...
        maxCdfError(j));
end

fprintf('\n');

%% ========================================================================
%  4. Final RQMC statistical consistency check
% ========================================================================

finalError = ...
    abs(terminalMean(end) - probCircle.F(end));

finalSE = terminalSE(end);

fprintf('Final terminal comparison:\n');
fprintf('  RQMC      = %.12e\n', terminalMean(end));
fprintf('  reference = %.12e\n', probCircle.F(end));
fprintf('  error     = %.12e\n', finalError);
fprintf('  SE        = %.12e\n', finalSE);

% Include a small deterministic floor because the replicate standard
% error itself is estimated from only J replicates.
acceptanceTolerance = max(5 * finalSE, 2e-5);

fprintf('  acceptance tolerance = %.12e\n\n', ...
    acceptanceTolerance);

assert(finalError < acceptanceTolerance, ...
    'Final RQMC terminal probability is inconsistent with reference.');

% Complete-CDF check using pointwise standard errors from the last run.
cdfTolerance = max( ...
    5 * max(rqmc.standardError), ...
    5e-5);

assert(maxCdfError(end) < cdfTolerance, ...
    'Final RQMC CDF is inconsistent with circular reference.');

%% ========================================================================
%  Same-rule local exceedance check
% ========================================================================

fprintf('Same-rule local exceedance check:\n');

localViolation = max( ...
    rqmc.pLocal - rqmc.F);

minimumMemoryCorrection = min( ...
    rqmc.memoryCorrection);

fprintf('  max[pLocal - F]       = %.12e\n', ...
    localViolation);

fprintf('  min memory correction = %.12e\n\n', ...
    minimumMemoryCorrection);

assert(localViolation < 1e-12, ...
    ['Local exceedance exceeds first-transition CDF despite using ', ...
     'the same RQMC directions.']);

assert(minimumMemoryCorrection > -1e-12, ...
    'The memory correction is materially negative.');

% At the inlet, local exceedance and first-transition probability must
% agree for every replicate.
inletDifference = max(abs( ...
    rqmc.pLocalReplicate(1, :) ...
    - rqmc.FReplicate(1, :)));

fprintf('  maximum inlet replicate difference = %.12e\n\n', ...
    inletDifference);

assert(inletDifference < 1e-12, ...
    'Local and first-transition probabilities differ at the inlet.');

%% ========================================================================
%  5. Reproducibility check
% ========================================================================

fprintf('Reproducibility check...\n');

rqmcRepeat = computeRQMCCDF( ...
    x, G, Gprime, r, ...
    angularCfg, ...
    cfg.maxDetection, ...
    cfg.numerics);

reproducibilityDefect = ...
    max(abs(rqmcRepeat.F - rqmc.F));

fprintf('  maximum repeated-run difference = %.12e\n\n', ...
    reproducibilityDefect);

assert(reproducibilityDefect == 0, ...
    'RQMC calculation is not reproducible under a fixed seed.');

%% ========================================================================
%  6. Replicate independence check
% ========================================================================

terminalSpread = std( ...
    rqmc.terminalProbabilityReplicate);

fprintf('Replicate spread:\n');
fprintf('  standard deviation among terminal estimates = %.12e\n\n', ...
    terminalSpread);

assert(terminalSpread > 0, ...
    'RQMC replicates appear to be identical.');

%% ========================================================================
%  7. Plots
% ========================================================================

figure('Color', 'w', ...
    'Name', 'RQMC convergence');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

errorbar( ...
    KValues, ...
    terminalMean, ...
    2 * terminalSE, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

hold on;

yline(probCircle.F(end), '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Circular reference');

set(gca, 'XScale', 'log');

xlabel('Directions per replicate K');
ylabel('F_{X_{tr}}(x_{max})');
title('RQMC terminal probability');
grid on;
box on;

nexttile;

loglog(KValues, maxCdfError, 'o-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'auto');

xlabel('Directions per replicate K');
ylabel('Maximum CDF error');
title('Error relative to circular reference');
grid on;
box on;

nexttile;
hold on;

plot(x, probCircle.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'Circular reference');

plot(x, rqmc.F, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'RQMC mean');

xlabel('x');
ylabel('F_{X_{tr}}(x)');
title('Complete CDF comparison');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, rqmc.standardError, ...
    'LineWidth', 1.5);

xlabel('x');
ylabel('Estimated standard error');
title('RQMC replicate uncertainty');
grid on;
box on;

sgtitle('Scrambled-Sobol angular integration');

fprintf('============================================================\n');
fprintf('All RQMC angular integration tests passed.\n');
fprintf('============================================================\n');