function result = computeProperComplexRQMCCDF( ...
    x, G, Gprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTEPROPERCOMPLEXRQMCCDF Proper-complex first-transition RQMC.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

J = angularCfg.numReplicates;

validateattributes(J, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, 'angularCfg.numReplicates');

FReplicate = zeros(Nx, J);
pLocalReplicate = zeros(Nx, J);

cfgLocal = struct();
cfgLocal.numerics = numericsCfg;
cfgLocal.maxDetection = maxDetectionCfg;
cfgLocal.gaussianConvention = 'properComplex';

replicateMeta = cell(J, 1);

for j = 1:J
    replicateCfg = angularCfg;
    replicateCfg.replicateIndex = j;

    [U, weights, angularMeta] = ...
        generateProperComplexAngularRule( ...
            r, replicateCfg);

    validateAngularRule(U, weights);

    switch lower(maxDetectionCfg.method)

        case 'grid'
            maxData = computeRunningMaxGrid( ...
                G, U, cfgLocal);

        case 'hermite'
            if isempty(Gprime)
                error('computeProperComplexRQMCCDF:MissingGprime', ...
                    'Hermite maximum detection requires Gprime.');
            end

            maxData = computeRunningMaxHermite( ...
                x, G, Gprime, U, cfgLocal);

        case 'adaptive'
            if isempty(Gprime)
                error('computeProperComplexRQMCCDF:MissingGprime', ...
                    'Adaptive maximum detection requires Gprime.');
            end

            maxData = computeRunningMaxAdaptive( ...
                x, G, Gprime, U, cfgLocal);

        otherwise
            error('computeProperComplexRQMCCDF:UnknownMaximumMethod', ...
                'Unknown maximum method "%s".', ...
                maxDetectionCfg.method);
    end

    prob = computeProperComplexTransitionCDF( ...
        maxData.m, weights, r);

    pLocal = ...
        computeProperComplexLocalExceedanceFromGains( ...
            maxData.a, weights, r);

    FReplicate(:, j) = prob.F;
    pLocalReplicate(:, j) = pLocal;

    replicateMeta{j} = struct();
    replicateMeta{j}.angularMeta = angularMeta;
    replicateMeta{j}.maximumInfo = maxData.info;
end

F = mean(FReplicate, 2);
S = 1 - F;

pLocal = mean(pLocalReplicate, 2);

memoryCorrectionReplicate = ...
    FReplicate - pLocalReplicate;

memoryCorrection = ...
    mean(memoryCorrectionReplicate, 2);

standardError = ...
    std(FReplicate, 0, 2) / sqrt(J);

pLocalStandardError = ...
    std(pLocalReplicate, 0, 2) / sqrt(J);

memoryCorrectionStandardError = ...
    std(memoryCorrectionReplicate, 0, 2) / sqrt(J);

pInterval = zeros(Nx, 1);
pInterval(1) = F(1);

if Nx > 1
    pInterval(2:end) = diff(F);
end

smallNegative = ...
    pInterval < 0 & pInterval > -1e-13;

pInterval(smallNegative) = 0;

localViolation = max(pLocal - F);

if localViolation > 1e-12
    error('computeProperComplexRQMCCDF:LocalViolation', ...
        ['The local probability exceeds the first-transition CDF ', ...
         'by %.3e.'], localViolation);
end

result = struct();

result.F = F;
result.S = S;
result.pLocal = pLocal;
result.memoryCorrection = memoryCorrection;

result.standardError = standardError;
result.pLocalStandardError = pLocalStandardError;
result.memoryCorrectionStandardError = ...
    memoryCorrectionStandardError;

result.FReplicate = FReplicate;
result.pLocalReplicate = pLocalReplicate;
result.memoryCorrectionReplicate = ...
    memoryCorrectionReplicate;

result.pInterval = pInterval;
result.pCensored = S(end);
result.terminalStandardError = standardError(end);

result.numDirections = ...
    angularCfg.numDirections;

result.numReplicates = J;
result.replicates = replicateMeta;

result.meta = struct();
result.meta.gaussianConvention = 'properComplex';
result.meta.complexDimension = r;

end
