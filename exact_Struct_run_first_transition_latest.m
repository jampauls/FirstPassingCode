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

clear; clc;

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
%   'matfile'       : load precomputed arrays from .mat file
%   'synthetic'     : use analytical/synthetic test problem
%   'ownsMarch'     : call OWNS propagator to generate B(x)
%   'ownsUserSelect': opens uigetfile to allow user to select file
cfg.dataSource = 'ownsmatfile';

% File path for precomputed data, if cfg.dataSource = 'matfile'.
%cfg.dataFile = fullfile('data', 'ninth_Run_i_10_j_1_solution.mat');
cfg.dataFile = fullfile('data', 'ninth_Run_i_1_j_1_solution.mat');

% -------------------------------------------------------------------------
% Interactive OWNS file selection
% -------------------------------------------------------------------------
if strcmp(cfg.dataSource, 'ownsUserSelect')
    % Initial directory shown by uigetfile when:
    %
    %   cfg.dataSource = 'ownsUserSelect'
    %
    % This may be an absolute or relative directory.
    cfg.owns.userSelectDirectory = ...
        '/data2/jampauls/DataforFigures/SaveData/FinalData';
    
    % Dialog title.
    cfg.owns.userSelectTitle = ...
        'Select an OWNS solution MAT-file';
    
    % Optional default filename. Leave empty to show all MAT-files.
    cfg.owns.userSelectDefaultFile = '';
    
    % If true, print the selected path.
    cfg.owns.userSelectVerbose = true;
end

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
cfg.angular.numDirections = 8192;

% RQMC options.
cfg.angular.numReplicates = 16;

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
% Transition-oriented stochastic dimension reduction
% -------------------------------------------------------------------------

cfg.reduction = struct();

% Options:
%   'none'
%   'givenBasis'
%   'integrated'
%   'eigUnion'
%   'hybrid'
%   'greedy'
cfg.reduction.method = 'integrated';

% Final reduced stochastic dimension.
cfg.reduction.targetRank = 30;

% Number of leading integrated-Gramian modes used to seed a hybrid basis.
cfg.reduction.numIntegratedModes = 10;

% Number of local eigenvectors collected at each selected station.
cfg.reduction.numModesPerStation = 3;

% Station-selection options for eigUnion/hybrid:
%   'records'    : stations where lambda_max(G) establishes a new record
%   'uniform'    : uniformly spaced stations
%   'all'        : every stored station
%   'recordsAndUniform'
cfg.reduction.stationSelection = 'recordsAndUniform';

% Number of uniformly selected stations, including endpoints.
cfg.reduction.numUniformStations = 12;

% Relative increase required to declare a new record station.
cfg.reduction.recordRelativeTolerance = 1e-6;

% Weight used in the integrated matrix:
%   'uniform'          : integral of G(x)
%   'traceNormalized'  : G(x)/trace(G(x))
%   'spectralNormalized': G(x)/lambda_max(G(x))
cfg.reduction.integrationWeight = 'uniform';

% Greedy enrichment options.
cfg.reduction.useGreedyEnrichment = true;
cfg.reduction.greedyTolerance = 1e-8;
cfg.reduction.maxGreedyModes = 50;

% Numerical rank tolerance when compressing a snapshot basis.
cfg.reduction.snapshotTolerance = 1e-12;

% Optional save of computed basis for further use
cfg.reduction.saveBasis = false;
if cfg.reduction.saveBasis
    cfg.reduction.saveName = 'ninth_Run_rank30_basis.mat';
end

% Optional user-supplied basis.
cfg.reduction.basisFile = '';

% Print diagnostics while constructing the basis.
cfg.reduction.verbose = true;

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
% OWNS data configuration
% -------------------------------------------------------------------------

cfg.owns = struct();

cfg.owns.dataFile = cfg.dataFile;%...
    %'/data2/jampauls/DataforFigures/SaveData/FirstPassingCode/data/ninth_Run_i_1_j_1_solution.mat';

cfg.owns.solutionVariable = 'solution';

% q is a synthesis factor:
%
%   qRealization(x_n) = squeeze(solution.q(:,n,:)) * w
cfg.owns.factorOrientation = 'synthesis';

% The first five variable blocks contribute to energy:
%
%   [rho, u, v, w, T]
%
% Pressure is the sixth block and is excluded.
cfg.owns.numEnergyVariables = 5;
cfg.owns.numStateVariables = 6;

cfg.owns.weightTolerance = 1e-12;

% Restrict the energy norm to a local delta99 region without modifying
% the propagated solution state. The smooth cosine taper avoids discrete
% energy jumps when the cutoff crosses a wall-normal grid row.
cfg.owns.trimToDelta99 = true;
cfg.owns.numDelta99 = 1;
cfg.owns.delta99MaskMethod = 'cosineTaper';% 'binary'
cfg.owns.delta99TaperWidth = 0.25;

% Print an inspection report while loading.
cfg.owns.verbose = true;

% -------------------------------------------------------------------------
% Streamwise-coordinate configuration
% -------------------------------------------------------------------------

cfg.coordinates = struct();

% Options:
%   'wallArcLength'
%   'referenceArcLength'
%   'wallX'
%   'referenceX'
%   'xi'
cfg.coordinates.method = 'wallArcLength';

% Wall-normal row used by reference methods.
% Row 1 is the wall for the current grid.
cfg.coordinates.referenceIndex = 1;

% Include z displacement in arc length, if available.
cfg.coordinates.includeZ = true;

% If true, arc length starts at zero.
cfg.coordinates.zeroOrigin = true;

% -------------------------------------------------------------------------
% Transition-threshold configuration
% -------------------------------------------------------------------------

cfg.threshold = struct();

% Options:
%   'auto'
%   'specifiedScalar'
%   'specifiedVector'
%   'initialMeanMultiple'
%   'localMeanMultiple'
%   'maxMeanFraction'
%   'initialQuantileApprox'
cfg.threshold.method = 'auto';

% Used by 'specifiedScalar'.
cfg.threshold.value = [];

% Used by 'specifiedVector'.
cfg.threshold.vector = [];

% Used by 'initialMeanMultiple' and 'localMeanMultiple'.
cfg.threshold.multiplier = 10;

% Used by 'maxMeanFraction'.
cfg.threshold.fraction = 1;

% Used by 'auto'. A constant threshold is set to this factor times the
% maximum mean energy over the OWNS domain.
cfg.threshold.autoMaxMeanFactor = 1;

% Used by 'initialQuantileApprox'. This is the desired inlet exceedance
% probability.
cfg.threshold.exceedanceProbability = 1e-3;

% Print a warning whenever a data-derived provisional threshold is used.
cfg.threshold.warnIfProvisional = true;

% -------------------------------------------------------------------------
% Threshold and inlet-amplitude exploration
% -------------------------------------------------------------------------

cfg.thresholdSweep = struct();

cfg.thresholdSweep.enable = false;

% Constant threshold values. If empty, values can be generated relative
% to the maximum mean energy.
cfg.thresholdSweep.values = [];

% Multiplicative factors applied to a reference threshold.
cfg.thresholdSweep.factors = ...
    logspace(-0.5, 0.5,  9);

% Inlet-factor amplitude multipliers:
%
%   B0(epsilon) = epsilon * B0(reference)
%
% Energy therefore scales by epsilon^2.
cfg.thresholdSweep.amplitudes = [0.5, 0.75, 1.0, 1.25, 1.5];

% If thresholdSweep.values is empty, use:
%
%   threshold = factor * max_x E[e(x)].
cfg.thresholdSweep.referenceMethod = 'maximumMeanEnergy';

% Optional threshold solution for a target terminal probability.
cfg.thresholdSweep.solveTarget.enable = false;
cfg.thresholdSweep.solveTarget.probability = 0.5;
cfg.thresholdSweep.solveTarget.amplitude = 1.0;
cfg.thresholdSweep.solveTarget.relativeTolerance = 1e-6;
cfg.thresholdSweep.solveTarget.probabilityTolerance = 1e-6;
cfg.thresholdSweep.solveTarget.maxIterations = 100;

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

B = [];
H = [];
A = [];
Vred = [];

switch lower(cfg.dataSource)

    case 'synthetic'
        fprintf('Generating synthetic test problem...\n');

        [x, B, H, eThresh, meta] = ...
            makeSyntheticProblem(cfg);

        B = convertFactorOrientation( ...
            B, cfg.factorOrientation);

        validateProblemData( ...
            x, B, H, eThresh, cfg);

        representation = 'factors';

        Nx = numel(x);
        [nq, rFull] = size(B{1});

    case 'matfile'
        fprintf('Loading generic problem data from: %s\n', ...
            cfg.dataFile);

        [x, B, H, eThresh, meta] = ...
            loadProblemData(cfg.dataFile, cfg);

        B = convertFactorOrientation( ...
            B, cfg.factorOrientation);

        validateProblemData( ...
            x, B, H, eThresh, cfg);

        representation = 'factors';

        Nx = numel(x);
        [nq, rFull] = size(B{1});

    case {'ownsmatfile', 'ownsuserselect'}
        
        if strcmpi(cfg.dataSource, 'ownsUserSelect')
            % Prompt the user for an OWNS solution file.
            selectedDataFile = ...
                selectOWNSDataFile(cfg.owns);
            
            % Store the selected path in cfg so it is preserved in saved
            % results and available to later code.
            cfg.owns.dataFile = selectedDataFile;
        else
            selectedDataFile = cfg.owns.dataFile;
        end
        
        problem = loadOWNSProblem( ...
            selectedDataFile, cfg);
        
        x = problem.x;
        A = problem.A;
        eThresh = problem.eThresh;
        meta = problem.meta;
        meta.selectedDataFile = selectedDataFile;
        
        representation = 'energyMatrices';
        
        Nx = numel(x);
        nq = meta.nq;
        rFull = size(A{1}, 1);
        
        clear problem;


    case 'ownsmarch'
        error('The direct OWNS marching interface is not yet implemented.');

    otherwise
        error('Unknown cfg.dataSource: %s', ...
            cfg.dataSource);
end

fprintf('Problem dimensions:\n');
fprintf('  Nx   = %d streamwise stations\n', Nx);
fprintf('  nq   = %d physical state dimension\n', nq);
fprintf('  r    = %d stochastic dimension before reduction\n', rFull);

%% ========================================================================
%  3. Optional transition-oriented dimension reduction
% ========================================================================

reductionInfo = struct();

switch lower(cfg.reduction.method)

    case 'none'
        fprintf('No stochastic dimension reduction.\n');

        r = rFull;

        if strcmpi(representation, 'factors')
            Buse = B;
            Ause = [];
        else
            Ause = A;
            Buse = [];
        end

    case 'givenbasis'
        fprintf('Loading user-supplied reduction basis...\n');

        tmp = load(cfg.reduction.basisFile);

        if ~isfield(tmp, 'Vred')
            error('run_first_transition:MissingVred', ...
                'The basis file must contain Vred.');
        end

        Vred = tmp.Vred;

        validateReductionBasis(Vred, rFull);

        if strcmpi(representation, 'factors')
            Buse = applyReductionToFactors(B, Vred);
            Ause = [];
        else
            Ause = applyReductionToMatrixFamily(A, Vred);
            Buse = [];
        end

        r = size(Vred, 2);

        reductionInfo.method = 'givenBasis';
        reductionInfo.targetRank = r;

    case {'integrated', 'eigunion', 'hybrid', 'greedy'}
        fprintf('Constructing transition-oriented basis...\n');

        if strcmpi(representation, 'factors')
            Gfull = buildGFamily( ...
                B, H, eThresh, cfg);
        else
            Gfull = buildGFromAFamily( ...
                A, eThresh, cfg);
        end

        [Vred, reductionInfo] = buildTransitionBasis( ...
            Gfull, x, cfg);

        validateReductionBasis(Vred, rFull);

        if strcmpi(representation, 'factors')
            Buse = applyReductionToFactors(B, Vred);
            Ause = [];
        else
            Ause = applyReductionToMatrixFamily(A, Vred);
            Buse = [];
        end

        r = size(Vred, 2);

        clear Gfull;

    otherwise
        error('Unknown reduction method: %s', ...
            cfg.reduction.method);
end

if cfg.reduction.saveBasis
    basisData = struct();

    basisData.Vred = Vred;
    basisData.reductionInfo = reductionInfo;
    basisData.thresholdMethod = cfg.threshold.method;
    basisData.threshold = eThresh;
    basisData.omega = meta.omega;
    basisData.beta = meta.beta;
    basisData.sourceFile = cfg.owns.dataFile;
    
    save( ...
        fullfile('data', cfg.reduction.saveName), ...
        '-struct', ...
        'basisData');
end

fprintf('Probability calculation stochastic dimension r = %d\n', r);

if ~strcmpi(cfg.reduction.method, 'none')
    fprintf('Reduction ratio: %d / %d = %.3f\n', ...
        r, rFull, r / rFull);
end

%% ========================================================================
%  4. Build normalized quadratic-form matrices G(x)
% ========================================================================

fprintf('Building normalized energy matrices G(x)...\n');

if strcmpi(representation, 'factors')
    G = buildGFamily( ...
        Buse, H, eThresh, cfg);
else
    G = buildGFromAFamily( ...
        Ause, eThresh, cfg);
end

validateGFamily(G, cfg);

Gprime = [];

if ismember(lower(cfg.maxDetection.method), ...
        {'hermite', 'adaptive'})

    fprintf('Building Gprime(x)...\n');

    Gprime = buildGprimeFamily( ...
        x, G, meta, cfg, Vred);

    validateGprimeFamily( ...
        x, G, Gprime);
end


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
    results.reductionInfo = reductionInfo;

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