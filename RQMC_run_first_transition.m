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
%  5--7. Angular integration, running maxima, and transition CDF
% ========================================================================

fprintf('Computing angular first-transition integral...\n');

if strcmpi(cfg.angular.method, 'rqmc') && ...
        cfg.angular.numReplicates > 1

    if cfg.angular.autoConverge
        prob = computeRQMCCDFConverged( ...
            x, G, Gprime, r, ...
            cfg.angular, ...
            cfg.maxDetection, ...
            cfg.numerics);
    else
        prob = computeRQMCCDF( ...
            x, G, Gprime, r, ...
            cfg.angular, ...
            cfg.maxDetection, ...
            cfg.numerics);
    end

    % No single angular rule represents the replicate-mean calculation.
    U = [];
    wAng = [];
    angularMeta = prob.meta;
    maxData = struct();
    maxData.info = struct( ...
        'method', cfg.maxDetection.method, ...
        'angularMethod', 'rqmcReplicateMean');

    fprintf('RQMC rule: K = %d directions per replicate, J = %d\n', ...
        prob.numDirections, prob.numReplicates);

    fprintf('Terminal RQMC standard error = %.6e\n', ...
        prob.terminalStandardError);

else
    [U, wAng, angularMeta] = ...
        generateAngularRule(r, cfg.angular);

    validateAngularRule(U, wAng);

    K = size(U, 2);
    fprintf('Angular rule: K = %d directions\n', K);

    fprintf('Computing directional running maxima...\n');

    switch lower(cfg.maxDetection.method)

        case 'grid'
            maxData = computeRunningMaxGrid(G, U, cfg);

        case 'hermite'
            maxData = computeRunningMaxHermite( ...
                x, G, Gprime, U, cfg);

        case 'adaptive'
            maxData = computeRunningMaxAdaptive( ...
                x, G, Gprime, U, cfg);

        otherwise
            error('Unknown max detection method: %s', ...
                cfg.maxDetection.method);
    end

    validateRunningMax(maxData.m);

    fprintf('Computing survival and transition CDF...\n');

    prob = computeTransitionCDF( ...
        maxData.m, wAng, r, cfg);
end

validateProbabilityCurves(prob, cfg);

fprintf('Terminal transition probability F(xmax) = %.6e\n', ...
    prob.F(end));

fprintf('Right-censoring probability S(xmax)     = %.6e\n', ...
    prob.S(end));

% A replicate-mean RQMC result does not have one primary angular rule.
% These variables will hold a separate rule used only for directional
% diagnostics, not for computing the primary probability curves.
Udiagnostic = [];
wDiagnostic = [];
diagnosticMaxData = [];

%% ========================================================================
%  8. Optional local exceedance probability check
% ========================================================================

if cfg.plot.showLocalProbability
    fprintf('Computing local exceedance probabilities for consistency...\n');

    if isfield(prob, 'pLocal')
        % For a replicate-mean RQMC calculation, computeRQMCCDF has
        % already evaluated the local exceedance probabilities using
        % the same angular directions as the first-transition CDF.
        pLocal = prob.pLocal;

        localCheckTol = 1e-12;

    else
        % A single deterministic or random angular rule is available.
        pLocal = computeLocalExceedance( ...
            G, U, wAng, r, cfg);

        localCheckTol = 1e-12;
    end

    maximumLocalViolation = max( ...
        pLocal(:) - prob.F(:));

    fprintf('Maximum pLocal - F_Xtr = %.3e\n', ...
        maximumLocalViolation);

    if maximumLocalViolation > localCheckTol
        warning('run_first_transition:LocalProbabilityCheck', ...
            ['Some local probabilities exceed the computed CDF by ', ...
             'more than the numerical tolerance.']);
    end
else
    pLocal = [];
end

%% ========================================================================
%  9. Post-processing: interval probabilities, quantiles, diagnostics
% ========================================================================

fprintf('Post-processing transition distribution...\n');

post = struct();

% -------------------------------------------------------------------------
% 9.1 Transition quantiles
% -------------------------------------------------------------------------

post.quantileLevels = ...
    [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95];

post.quantiles = computeTransitionQuantiles( ...
    x, prob.F, post.quantileLevels);

% -------------------------------------------------------------------------
% 9.2 Memory correction and its optional RQMC uncertainty
% -------------------------------------------------------------------------

post.memoryCorrection = [];
post.memoryCorrectionStandardError = [];

if isfield(prob, 'memoryCorrection')
    % RQMC replicate-paired estimate.
    post.memoryCorrection = ...
        prob.memoryCorrection;

    if isfield(prob, 'memoryCorrectionStandardError')
        post.memoryCorrectionStandardError = ...
            prob.memoryCorrectionStandardError;
    end

elseif ~isempty(pLocal)
    % Single-rule deterministic or random angular calculation.
    post.memoryCorrection = ...
        prob.F(:) - pLocal(:);
end

% -------------------------------------------------------------------------
% 9.3 Directional diagnostics
% -------------------------------------------------------------------------

if ~isempty(U)
    % A single primary angular rule was used, so reuse it for diagnostics.
    post.directional = computeDirectionalDiagnostics( ...
        maxData.m, ...
        U, ...
        wAng, ...
        r, ...
        cfg);

else
    % The primary result is an average over multiple RQMC replicates and
    % therefore does not have one direction matrix U. Generate a separate
    % angular rule used only for directional diagnostics.
    diagnosticAngularCfg = cfg.angular;
    diagnosticAngularCfg.method = 'rqmc';

    % Use a scrambling not included in the primary replicate mean.
    diagnosticAngularCfg.replicateIndex = ...
        cfg.angular.numReplicates + 1;

    [Udiagnostic, wDiagnostic, diagnosticAngularMeta] = ...
        generateAngularRule( ...
            r, diagnosticAngularCfg);

    validateAngularRule( ...
        Udiagnostic, wDiagnostic);

    switch lower(cfg.maxDetection.method)

        case 'grid'
            diagnosticMaxData = ...
                computeRunningMaxGrid( ...
                    G, Udiagnostic, cfg);

        case 'hermite'
            if isempty(Gprime)
                error('run_first_transition:MissingGprime', ...
                    ['Hermite directional diagnostics require ', ...
                     'Gprime.']);
            end

            diagnosticMaxData = ...
                computeRunningMaxHermite( ...
                    x, G, Gprime, ...
                    Udiagnostic, cfg);

        case 'adaptive'
            if isempty(Gprime)
                error('run_first_transition:MissingGprime', ...
                    ['Adaptive directional diagnostics require ', ...
                     'Gprime.']);
            end

            diagnosticMaxData = ...
                computeRunningMaxAdaptive( ...
                    x, G, Gprime, ...
                    Udiagnostic, cfg);

        otherwise
            error('run_first_transition:UnknownMaximumMethod', ...
                'Unknown maximum method "%s".', ...
                cfg.maxDetection.method);
    end

    validateRunningMax( ...
        diagnosticMaxData.m);

    post.directional = computeDirectionalDiagnostics( ...
        diagnosticMaxData.m, ...
        Udiagnostic, ...
        wDiagnostic, ...
        r, ...
        cfg);

    post.directional.angularMeta = ...
        diagnosticAngularMeta;

    post.directional.isPrimaryProbabilityRule = false;
end


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