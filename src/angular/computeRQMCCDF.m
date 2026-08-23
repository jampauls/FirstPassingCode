function rqmc = computeRQMCCDF( ...
    x, G, Gprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTERQMCCDF Evaluate first-transition probabilities with RQMC.
%
% The first-transition CDF and local exceedance probability are evaluated
% using the same angular directions in every replicate. This preserves
%
%   pLocal(x_n) <= F_Xtr(x_n)
%
% replicate by replicate, up to floating-point error.
%
% MATLAB version: R2020b

if ~strcmpi(angularCfg.method, 'rqmc')
    error('computeRQMCCDF:InvalidMethod', ...
        'angularCfg.method must be ''rqmc''.');
end

J = angularCfg.numReplicates;

validateattributes(J, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, 'angularCfg.numReplicates');

Nx = numel(x);

FReplicate = zeros(Nx, J);
pLocalReplicate = zeros(Nx, J);
terminalProbability = zeros(J, 1);
replicateMeta = cell(J, 1);

maximumDirectionalViolation = zeros(J, 1);

for j = 1:J
    replicateCfg = angularCfg;
    replicateCfg.replicateIndex = j;

    replicate = evaluateAngularReplicate( ...
        x, G, Gprime, r, ...
        replicateCfg, maxDetectionCfg, numericsCfg);

    FReplicate(:, j) = replicate.prob.F;
    terminalProbability(j) = replicate.prob.F(end);

    % The endpoint directional gains and running maxima are computed
    % using exactly the same angular directions.
    if isfield(replicate.maxData, 'a') && ...
            ~isempty(replicate.maxData.a)

        pLocalReplicate(:, j) = ...
            computeLocalExceedanceFromGains( ...
                replicate.maxData.a, ...
                replicate.weights, ...
                r);

        % Check the stronger direction-by-direction condition:
        %
        %   a(k,n) <= m(k,n).
        maximumDirectionalViolation(j) = max( ...
            replicate.maxData.a(:) ...
            - replicate.maxData.m(:));

    else
        % Fallback for a future maximum implementation that does not
        % retain endpoint gains.
        pLocalReplicate(:, j) = ...
            computeLocalExceedance( ...
                G, ...
                replicate.U, ...
                replicate.weights, ...
                r, ...
                struct());

        maximumDirectionalViolation(j) = NaN;
    end

    replicateMeta{j} = struct( ...
        'angularMeta', replicate.angularMeta, ...
        'maxDataInfo', replicate.maxData.info);
end

% -------------------------------------------------------------------------
% Replicate means
% -------------------------------------------------------------------------

F = mean(FReplicate, 2);
S = 1 - F;

pLocal = mean(pLocalReplicate, 2);

% -------------------------------------------------------------------------
% Replicate uncertainty estimates
% -------------------------------------------------------------------------

standardDeviation = std(FReplicate, 0, 2);
standardError = standardDeviation / sqrt(J);

pLocalStandardDeviation = std(pLocalReplicate, 0, 2);
pLocalStandardError = ...
    pLocalStandardDeviation / sqrt(J);

% The replicate-level differences are useful because the two estimators
% are correlated when computed with the same directions.
memoryCorrectionReplicate = ...
    FReplicate - pLocalReplicate;

memoryCorrection = ...
    mean(memoryCorrectionReplicate, 2);

memoryCorrectionStandardDeviation = ...
    std(memoryCorrectionReplicate, 0, 2);

memoryCorrectionStandardError = ...
    memoryCorrectionStandardDeviation / sqrt(J);

% -------------------------------------------------------------------------
% Interval transition probabilities
% -------------------------------------------------------------------------

pInterval = zeros(Nx, 1);
pInterval(1) = F(1);

if Nx > 1
    pInterval(2:end) = diff(F);
end

smallNegative = ...
    pInterval < 0 & pInterval > -1e-13;

pInterval(smallNegative) = 0;

if any(pInterval < -1e-12)
    error('computeRQMCCDF:NonmonotoneMeanCDF', ...
        'The replicate-mean CDF is materially nonmonotone.');
end

pCensored = S(end);

% -------------------------------------------------------------------------
% Discrete hazard
% -------------------------------------------------------------------------

hDiscrete = nan(Nx, 1);
hDiscrete(1) = F(1);

for n = 2:Nx
    if S(n - 1) > 0
        hDiscrete(n) = ...
            (S(n - 1) - S(n)) / S(n - 1);
    end
end

% -------------------------------------------------------------------------
% Consistency checks
% -------------------------------------------------------------------------

localViolation = max(pLocal - F);

if localViolation > 1e-12
    error('computeRQMCCDF:LocalProbabilityViolation', ...
        ['The same-rule local probability exceeds the first-transition ', ...
         'CDF by %.3e.'], localViolation);
end

finiteDirectionalViolation = ...
    maximumDirectionalViolation( ...
        isfinite(maximumDirectionalViolation));

if ~isempty(finiteDirectionalViolation)
    maxDirectionalViolation = ...
        max(finiteDirectionalViolation);

    if maxDirectionalViolation > 1e-12
        error('computeRQMCCDF:RunningMaximumViolation', ...
            ['An endpoint directional gain exceeds its running maximum ', ...
             'by %.3e.'], maxDirectionalViolation);
    end
else
    maxDirectionalViolation = NaN;
end

% -------------------------------------------------------------------------
% Output structure
% -------------------------------------------------------------------------

rqmc = struct();

% Core probability fields.
rqmc.F = F;
rqmc.S = S;
rqmc.pInterval = pInterval;
rqmc.pCensored = pCensored;
rqmc.hDiscrete = hDiscrete;

% Local exceedance and memory correction.
rqmc.pLocal = pLocal;
rqmc.memoryCorrection = memoryCorrection;

% RQMC uncertainty output.
rqmc.standardDeviation = standardDeviation;
rqmc.standardError = standardError;

rqmc.pLocalStandardDeviation = ...
    pLocalStandardDeviation;

rqmc.pLocalStandardError = ...
    pLocalStandardError;

rqmc.memoryCorrectionStandardDeviation = ...
    memoryCorrectionStandardDeviation;

rqmc.memoryCorrectionStandardError = ...
    memoryCorrectionStandardError;

% Replicate data.
rqmc.FReplicate = FReplicate;
rqmc.pLocalReplicate = pLocalReplicate;
rqmc.memoryCorrectionReplicate = ...
    memoryCorrectionReplicate;

rqmc.terminalProbabilityReplicate = ...
    terminalProbability;

rqmc.terminalStandardError = ...
    standardError(end);

rqmc.numDirections = ...
    angularCfg.numDirections;

rqmc.numReplicates = J;
rqmc.replicates = replicateMeta;

rqmc.meta = struct();
rqmc.meta.method = 'rqmc';
rqmc.meta.totalProbability = ...
    sum(pInterval) + pCensored;

rqmc.meta.maximumLocalCdfViolation = ...
    localViolation;

rqmc.meta.maximumDirectionalViolation = ...
    maxDirectionalViolation;

end