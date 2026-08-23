% run_first_transition.m
%
% Main driver for memory-aware first-transition prediction from a
% propagated stochastic factor B(x).
%
% Core model:
%   q(x) = B(x) w,       w ~ N(0,I_r)     [default real convention]
%   e(x) = q(x)^* H(x) q(x)
%   G(x) = real(B(x)^* H(x) B(x)) / e_thres(x)
%
% First-transition CDF:
%   F_Xtr(x) = E_u[ chi2sf( 1 / max_{xi<=x} u^T G(xi) u ; r ) ]
%
% MATLAB version target: R2020b

clear; %clc;

%% ========================================================================
%  0. Path setup
% ========================================================================

% Add project subdirectories as they are created.
% Keep this section minimal at first.
addpath(genpath(fullfile(pwd, 'src')));
addpath(genpath(fullfile(pwd, 'examples')));
addpath(genpath(fullfile(pwd, 'data')));
addpath(genpath(fullfile(pwd, 'tests')));
addpath(genpath(fullfile(pwd, 'figures')));

%% ========================================================================
%  1. User configuration
% ========================================================================

cfg = struct();

% -------------------------------------------------------------------------
% 1.1 Problem convention
% -------------------------------------------------------------------------

% Gaussian convention:
%   'real'          : w ~ N(0,I_r), radial variable R^2 ~ chi2(r)
%   'properComplex' : z ~ CN(0,I_r), radial variable |z|^2 ~ Gamma(r,1)
%
% Start with 'real'. Proper-complex support can be added later.
cfg.gaussianConvention = 'real';

% Factor orientation:
%   'synthesis' : q = B*w, C = B*B'
%   'row'       : C = B_row'*B_row, convert internally via B = B_row'
cfg.factorOrientation = 'synthesis';

% Whether the physical OWNS/Fourier state is complex.
cfg.isPhysicalStateComplex = true;

% -------------------------------------------------------------------------
% 1.2 Data source
% -------------------------------------------------------------------------

% Select how B(x), H(x), threshold, and x-grid will be obtained.
%
% Options to implement later:
%   'matfile'      : load precomputed arrays from .mat file
%   'synthetic'    : use analytical/synthetic test problem
%   'ownsMarch'    : call OWNS propagator to generate B(x)
cfg.dataSource = 'synthetic';

% File path for precomputed data, if cfg.dataSource = 'matfile'.
cfg.dataFile = fullfile('data', 'example_B_H_threshold.mat');

% -------------------------------------------------------------------------
% 1.3 Streamwise maximum detection
% -------------------------------------------------------------------------

% Maximum detection method:
%   'grid'       : endpoint-only maxima over stored x stations
%   'hermite'    : cubic Hermite interval maxima using G and G'
%   'adaptive'   : adaptive refinement, to be implemented later
cfg.maxDetection.method = 'grid';

% Tolerances for later adaptive maximum detection.
cfg.maxDetection.absTol = 1e-8;
cfg.maxDetection.relTol = 1e-6;
cfg.maxDetection.maxRefineLevel = 10;
cfg.maxDetection.method = 'hermite';

% -------------------------------------------------------------------------
% 1.4 Angular integration
% -------------------------------------------------------------------------

% Angular integration method:
%   'random'     : iid random directions on sphere
%   'rqmc'       : Sobol/RQMC directions, later
%   'circle'     : deterministic rule for r = 2
%   'user'       : supplied directions and weights

cfg.angular.method = 'rqmc';

% cfg.angular.numDirections = 5000;
cfg.angular.numDirections = 4096;

% RQMC options.
cfg.angular.numReplicates = 8;

% Independent scrambling seed for reproducibility.
cfg.angular.randomSeed = 1;

% Sobol' construction options. Use integer values in R2020b.
cfg.angular.sobolSkip = 1024;
cfg.angular.sobolLeap = 0;

% Whether to apply Matousek--Affine--Owen scrambling.
cfg.angular.scramble = true;

% Avoid evaluating the inverse normal CDF at exactly zero or one.
cfg.angular.probabilityClip = 1e-12;

% Optional automatic angular-convergence study.
cfg.angular.autoConverge = false;
cfg.angular.initialNumDirections = 256;
cfg.angular.maxNumDirections = 16384;
cfg.angular.cdfTolerance = 1e-5;
cfg.angular.standardErrorTolerance = 1e-5;

% -------------------------------------------------------------------------
% 1.5 Dimension reduction
% -------------------------------------------------------------------------

% Reduction mode:
%   'none'       : use full stochastic dimension r
%   'givenBasis' : use user-supplied V_s
%   'eigUnion'   : construct transition-oriented basis, later
%   'greedy'     : greedy transition-oriented basis, later
cfg.reduction.method = 'none';
cfg.reduction.targetRank = [];
cfg.reduction.basisFile = '';

% -------------------------------------------------------------------------
% 1.6 Numerical options
% -------------------------------------------------------------------------

cfg.numerics.symmetrizeG = true;
cfg.numerics.clipSmallNegativeGains = true;
cfg.numerics.psdProjection = false;
cfg.numerics.psdTol = 1e-10;

% Chi-square tail evaluation:
% MATLAB has chi2cdf in Statistics Toolbox.
% Upper tail can be evaluated with chi2cdf(z,r,'upper') in many MATLAB
% versions. Confirm compatibility in R2020b.
cfg.numerics.useUpperTail = true;

% -------------------------------------------------------------------------
% 1.7 Output and plotting
% -------------------------------------------------------------------------

cfg.output.saveResults = false;
cfg.output.resultsFile = fullfile('data', 'results_first_transition.mat');

cfg.plot.enable = true;
cfg.plot.showLocalProbability = true;
cfg.plot.showIntervalProbabilities = true;

%% ========================================================================
%  2. Load or generate problem data
% ========================================================================

% Required data after this section:
%
%   x           : [Nx x 1] streamwise stations
%   B           : propagated synthesis factors
%                 recommended storage: cell array B{n}, each [nq x r]
%   H           : energy matrices
%                 recommended storage: cell array H{n}, each [nq x nq]
%   eThresh     : [Nx x 1] threshold values
%
% Optional:
%   Gprime      : cell array of G'(x_n), for Hermite/adaptive maxima
%   meta        : struct with information about dimensions, units, etc.

switch lower(cfg.dataSource)

    case 'synthetic'
        fprintf('Generating synthetic test problem...\n');
        %cfg.synthetic.Nx = 41;
        [x, B, H, eThresh, meta] = makeSyntheticProblem(cfg);

    case 'matfile'
        fprintf('Loading problem data from: %s\n', cfg.dataFile);
        [x, B, H, eThresh, meta] = loadProblemData(cfg.dataFile, cfg);

    case 'ownsmarch'
        fprintf('Calling OWNS propagation workflow...\n');
        [x, B, H, eThresh, meta] = runOWNSWorkflow(cfg);

    otherwise
        error('Unknown cfg.dataSource: %s', cfg.dataSource);
end

% Convert all factors to the synthesis convention used internally:
%   q(x) = B(x)w,    C(x,x) = B(x)B(x)^*
B = convertFactorOrientation(B, cfg.factorOrientation);

% Basic dimension checks.
validateProblemData(x, B, H, eThresh, cfg);

Nx = numel(x);
[nq, rFull] = size(B{1});

fprintf('Problem dimensions:\n');
fprintf('  Nx   = %d streamwise stations\n', Nx);
fprintf('  nq   = %d physical state dimension\n', nq);
fprintf('  r    = %d stochastic dimension before reduction\n', rFull);

%% ========================================================================
%  3. Optional transition-oriented dimension reduction
% ========================================================================

% After this section:
%   Buse : cell array of propagated factors to use
%   r    : stochastic dimension actually used in probability calculation
%   Vred : reduction basis, if applicable

switch lower(cfg.reduction.method)

    case 'none'
        fprintf('No stochastic dimension reduction.\n');
        Buse = B;
        Vred = [];
        r = rFull;

    case 'givenbasis'
        fprintf('Loading user-supplied reduction basis...\n');
        tmp = load(cfg.reduction.basisFile);
        Vred = tmp.Vred;
        Buse = applyReductionToFactors(B, Vred);
        r = size(Vred, 2);

    case {'eigunion', 'greedy'}
        fprintf('Constructing transition-oriented basis...\n');

        % For these methods, we typically need preliminary G matrices.
        Gfull = buildGFamily(B, H, eThresh, cfg);

        Vred = buildTransitionBasis(Gfull, x, cfg);
        Buse = applyReductionToFactors(B, Vred);
        r = size(Vred, 2);

    otherwise
        error('Unknown cfg.reduction.method: %s', cfg.reduction.method);
end

fprintf('Probability calculation stochastic dimension r = %d\n', r);

%% ========================================================================
%  4. Build normalized quadratic-form matrices G(x)
% ========================================================================

fprintf('Building normalized energy matrices G(x)...\n');

G = buildGFamily(Buse, H, eThresh, cfg);

% Optional derivative matrices, not implemented initially.
Gprime = [];

if ismember(lower(cfg.maxDetection.method), {'hermite', 'adaptive'})
    fprintf('Building Gprime(x)...\n');

    Gprime = buildGprimeFamily( ...
        x, G, meta, cfg, Vred);
end

validateGFamily(G, cfg);
validateGprimeFamily(x, G, Gprime);

%% ========================================================================
%  5. Generate angular directions and weights
% ========================================================================

fprintf('Generating angular directions...\n');

[U, wAng, angularMeta] = generateAngularRule(r, cfg.angular);

% U storage convention:
%   U is [r x K], columns are unit directions.
%
% wAng:
%   [K x 1], nonnegative, sum to 1.

validateAngularRule(U, wAng);

K = size(U, 2);
fprintf('Angular rule: K = %d directions\n', K);

%% ========================================================================
%  6. Compute directional running maxima
% ========================================================================

fprintf('Computing directional running maxima...\n');

switch lower(cfg.maxDetection.method)

    case 'grid'
        % Endpoint-only streamwise maxima.
        maxData = computeRunningMaxGrid(G, U, cfg);

    case 'hermite'
        % Cubic Hermite maxima on intervals using G and Gprime.
        maxData = computeRunningMaxHermite(x, G, Gprime, U, cfg);

    case 'adaptive'
        % Adaptive streamwise refinement, to be added later.
        maxData = computeRunningMaxAdaptive(x, G, Gprime, U, cfg);

    otherwise
        error('Unknown max detection method: %s', cfg.maxDetection.method);
end

% Expected fields in maxData:
%   maxData.m        : [K x Nx] running maxima m_k(x_n)
%   maxData.a        : optional [K x Nx] endpoint gains
%   maxData.info     : diagnostics

validateRunningMax(maxData.m);

%% ========================================================================
%  7. Compute survival function and first-transition CDF
% ========================================================================

fprintf('Computing survival and transition CDF...\n');

prob = computeTransitionCDF(maxData.m, wAng, r, cfg);

% Expected fields in prob:
%   prob.S           : [Nx x 1] survival probability
%   prob.F           : [Nx x 1] transition CDF
%   prob.pInterval   : [Nx x 1] p0 and interval probabilities
%   prob.pCensored   : scalar
%   prob.hDiscrete   : [Nx x 1] optional discrete hazard
%   prob.meta        : diagnostics

validateProbabilityCurves(prob, cfg);

fprintf('Terminal transition probability F(xmax) = %.6e\n', prob.F(end));
fprintf('Right-censoring probability S(xmax)     = %.6e\n', prob.S(end));

%% ========================================================================
%  8. Optional local exceedance probability check
% ========================================================================

if cfg.plot.showLocalProbability
    fprintf('Computing local exceedance probabilities for consistency...\n');

    pLocal = computeLocalExceedance(G, U, wAng, r, cfg);

    % Check pLocal(x) <= F_Xtr(x), up to tolerance.
    localCheckTol = 1e-8;
    if any(pLocal(:) > prob.F(:) + localCheckTol)
        warning('Some local exceedance probabilities exceed first-transition CDF.');
    end
else
    pLocal = [];
end

%% ========================================================================
%  9. Post-processing: interval probabilities, quantiles, diagnostics
% ========================================================================

fprintf('Post-processing transition distribution...\n');

post = struct();

post.quantileLevels = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
post.quantiles = computeTransitionQuantiles(x, prob.F, post.quantileLevels);

post.memoryCorrection = [];
if ~isempty(pLocal)
    post.memoryCorrection = prob.F(:) - pLocal(:);
end

% Directional contribution diagnostics.
post.directional = computeDirectionalDiagnostics(maxData.m, U, wAng, r, cfg);

%% ========================================================================
%  10. Plot results
% ========================================================================

if cfg.plot.enable
    fprintf('Plotting results...\n');

    plotTransitionResults(x, prob, pLocal, post, cfg);
end

if strcmpi(cfg.dataSource, 'synthetic') && cfg.plot.enable
    plotSyntheticProblem(x, meta);
end

%% ========================================================================
%  11. Save results
% ========================================================================

if cfg.output.saveResults
    fprintf('Saving results to: %s\n', cfg.output.resultsFile);

    results = struct();
    results.cfg = cfg;
    results.x = x;
    results.prob = prob;
    results.post = post;
    results.meta = meta;
    results.angularMeta = angularMeta;
    results.maxDataInfo = maxData.info;

    if ~isempty(Vred)
        results.Vred = Vred;
    end

    save(cfg.output.resultsFile, 'results', '-v7.3');
end

fprintf('Done.\n');

%% ========================================================================
%  Local notes
% ========================================================================
%
% Functions to implement first:
%
%   makeSyntheticProblem.m
%   validateProblemData.m
%   buildGFamily.m
%   generateAngularRule.m
%   computeRunningMaxGrid.m
%   computeTransitionCDF.m
%   computeLocalExceedance.m
%   computeTransitionQuantiles.m
%   plotTransitionResults.m
%
% Functions to implement later:
%
%   loadProblemData.m
%   runOWNSWorkflow.m
%   buildGprimeFamily.m
%   computeRunningMaxHermite.m
%   computeRunningMaxAdaptive.m
%   buildTransitionBasis.m
%   applyReductionToFactors.m
%   computeDirectionalDiagnostics.m