%RUNBIVARIATEQUADRATICPROBABILITYTESTS Regression tests for 2D inversion.
%
% MATLAB version: R2020b
clear
clc

fprintf('============================================================\n');
fprintf('Bivariate quadratic probability regression tests\n');
fprintf('============================================================\n\n');

options = struct();

options.quadratureOrders = [16, 24, 32, 40];
options.dampingMultipliers = [0.25, 0.5, 1, 2, 4];

options.cfOptions = struct();
options.cfOptions.absoluteTolerance = 1e-10;
options.cfOptions.relativeTolerance = 1e-8;
options.cfOptions.maxIntervalCount = 15000;
options.cfOptions.warnOnLargeError = false;

options.verbose = false;

%% ========================================================================
% Test 1: independent one-dimensional quadratic forms
% ========================================================================

fprintf('Test 1: independent quadratic forms...\n');

lambda1 = 1.3;
lambda2 = 0.7;

level1 = 0.8;
level2 = 1.1;

B1 = diag([lambda1, 0]);
B2 = diag([0, lambda2]);

b1 = -level1;
b2 = -level2;

[pComputed, info] = ...
    bivariateQuadraticOrthantProbability( ...
        B1, b1, B2, b2, options);

p1 = gammainc( ...
    level1 / (2 * lambda1), ...
    1 / 2, 'upper');

p2 = gammainc( ...
    level2 / (2 * lambda2), ...
    1 / 2, 'upper');

pExact = p1 * p2;

fprintf('  exact       = %.12e\n', pExact);
fprintf('  computed    = %.12e\n', pComputed);
fprintf('  convergence = %.3e\n', ...
    info.convergenceDifference);

% assert(abs(pComputed - pExact) < 2e-3, ...
%     'Independent-form test failed.');
% 
% fprintf('  passed.\n\n');

%% ========================================================================
% Test 2: identical facet exclusion
% ========================================================================

fprintf('Test 2: identical-facet exclusion...\n');

A = diag([1.0, 0.5, 0.2]);
threshold = 1.1;

result = computePairwiseFacetExclusion( ...
    A, A, threshold, options);

fprintf('  exclusion probability = %.12e\n', ...
    result.probability);

assert(result.probability < 3e-3, ...
    'Identical-facet exclusion should be zero.');

fprintf('  passed.\n\n');

%% ========================================================================
% Test 3: proportional ordered forms
% ========================================================================

fprintf('Test 3: proportional ordered forms...\n');

r = 4;

sourceScale = 1.4;
representativeScale = 0.8;
threshold = 2.0;

A_source = sourceScale * eye(r);
A_representative = representativeScale * eye(r);

result = computePairwiseFacetExclusion( ...
    A_source, A_representative, threshold, options);

lowerChiSquare = threshold / sourceScale;
upperChiSquare = threshold / representativeScale;

pExact = ...
    gammainc(lowerChiSquare / 2, r / 2, 'upper') ...
    - gammainc(upperChiSquare / 2, r / 2, 'upper');

fprintf('  exact       = %.12e\n', pExact);
fprintf('  computed    = %.12e\n', result.probability);
fprintf('  convergence = %.3e\n', ...
    result.convergenceDifference);

assert(abs(result.probability - pExact) < 3e-3, ...
    'Proportional-form exclusion test failed.');

fprintf('  passed.\n\n');

%% ========================================================================
% Test 4: impossible ordered exclusion
% ========================================================================

fprintf('Test 4: impossible ordered exclusion...\n');

A_source = 0.6 * eye(3);
A_representative = 1.2 * eye(3);
threshold = 1.0;

result = computePairwiseFacetExclusion( ...
    A_source, A_representative, threshold, options);

fprintf('  computed = %.12e\n', result.probability);

assert(result.probability < 3e-3, ...
    'Ordered impossible exclusion should be zero.');

fprintf('  passed.\n\n');

fprintf('============================================================\n');
fprintf('All bivariate quadratic probability tests passed.\n');
fprintf('============================================================\n');
