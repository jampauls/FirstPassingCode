function study = computeCoupledRankStudy( ...
    x, Gfull, GprimeFull, Vnested, rankValues, ...
    angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTECOUPLEDRANKSTUDY Compare nested ranks with common RQMC coordinates.
%
% Each RQMC replicate is generated once in the full stochastic dimension.
% Rank s uses the first s coordinates and first s columns of Vnested.
%
% This produces correlated rank estimates and permits paired uncertainty
% estimates for CDF differences.
%
% Inputs:
%   x               Nx-by-1 streamwise coordinate
%   Gfull           cell array of rFull-by-rFull matrices
%   GprimeFull      derivative matrices, or []
%   Vnested         rFull-by-rFull ordered orthonormal basis
%   rankValues      tested nested ranks
%   angularCfg      RQMC configuration
%   maxDetectionCfg maximum-detection configuration
%   numericsCfg     numerical configuration
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);
rFull = size(Gfull{1}, 1);

validateReductionBasis(Vnested, rFull);

rankValues = unique(round(rankValues(:).'));

if any(rankValues < 1) || any(rankValues > rFull)
    error('computeCoupledRankStudy:InvalidRanks', ...
        'All tested ranks must lie between 1 and rFull.');
end

K = angularCfg.numDirections;
J = angularCfg.numReplicates;
numRanks = numel(rankValues);

validateattributes(J, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, 'angularCfg.numReplicates');

weights = ones(K, 1) / K;

% Dimensions:
%
%   FReplicate(:,rankIndex,replicateIndex)
FReplicate = zeros(Nx, numRanks, J);

pLocalReplicate = zeros(Nx, numRanks, J);

% Project matrix families once per tested rank.
Gred = cell(numRanks, 1);
GprimeRed = cell(numRanks, 1);
reductionDiagnostics = cell(numRanks, 1);

for jRank = 1:numRanks
    s = rankValues(jRank);
    Vs = Vnested(:, 1:s);

    Gred{jRank} = projectMatrixFamily(Gfull, Vs);

    if ~isempty(GprimeFull)
        GprimeRed{jRank} = ...
            projectMatrixFamily(GprimeFull, Vs);
    else
        GprimeRed{jRank} = [];
    end

    reductionDiagnostics{jRank} = ...
        computeReductionDiagnostics(Gfull, Vs);
end

fprintf('\nCoupled RQMC rank calculation:\n');
fprintf('  full dimension       = %d\n', rFull);
fprintf('  directions/replicate = %d\n', K);
fprintf('  replicates           = %d\n\n', J);

for jRep = 1:J
    fprintf('  replicate %d of %d\n', jRep, J);

    Zfull = generateScrambledSobolNormals( ...
        rFull, K, jRep, angularCfg);

    for jRank = 1:numRanks
        s = rankValues(jRank);

        U = nestedNormalsToDirections(Zfull, s);

        cfgLocal = struct();
        cfgLocal.numerics = numericsCfg;
        cfgLocal.maxDetection = maxDetectionCfg;
        cfgLocal.gaussianConvention = 'real';

        switch lower(maxDetectionCfg.method)

            case 'grid'
                maxData = computeRunningMaxGrid( ...
                    Gred{jRank}, U, cfgLocal);

            case 'hermite'
                if isempty(GprimeRed{jRank})
                    error('computeCoupledRankStudy:MissingGprime', ...
                        'Hermite maximum detection requires Gprime.');
                end

                maxData = computeRunningMaxHermite( ...
                    x, ...
                    Gred{jRank}, ...
                    GprimeRed{jRank}, ...
                    U, ...
                    cfgLocal);

            case 'adaptive'
                if isempty(GprimeRed{jRank})
                    error('computeCoupledRankStudy:MissingGprime', ...
                        'Adaptive maximum detection requires Gprime.');
                end

                maxData = computeRunningMaxAdaptive( ...
                    x, ...
                    Gred{jRank}, ...
                    GprimeRed{jRank}, ...
                    U, ...
                    cfgLocal);

            otherwise
                error('computeCoupledRankStudy:UnknownMaximumMethod', ...
                    'Unknown maximum method "%s".', ...
                    maxDetectionCfg.method);
        end

        prob = computeTransitionCDF( ...
            maxData.m, weights, s, cfgLocal);

        pLocal = computeLocalExceedanceFromGains( ...
            maxData.a, weights, s);

        FReplicate(:, jRank, jRep) = prob.F;
        pLocalReplicate(:, jRank, jRep) = pLocal;
    end
end

% -------------------------------------------------------------------------
% Rank-wise means and standard errors
% -------------------------------------------------------------------------

FMean = mean(FReplicate, 3);
FStandardDeviation = std(FReplicate, 0, 3);
FStandardError = FStandardDeviation / sqrt(J);

pLocalMean = mean(pLocalReplicate, 3);
pLocalStandardError = ...
    std(pLocalReplicate, 0, 3) / sqrt(J);

% -------------------------------------------------------------------------
% Paired differences between successive tested ranks
% -------------------------------------------------------------------------

successiveDifferenceMean = nan(Nx, numRanks);
successiveDifferenceSE = nan(Nx, numRanks);
maximumSuccessiveDifference = nan(numRanks, 1);
maximumSuccessiveDifferenceSE = nan(numRanks, 1);
maximumSuccessiveDifferenceIndex = nan(numRanks, 1);

for jRank = 2:numRanks
    differenceReplicate = ...
        squeeze(FReplicate(:, jRank, :) ...
        - FReplicate(:, jRank - 1, :));

    differenceMean = mean(differenceReplicate, 2);
    differenceSE = std(differenceReplicate, 0, 2) / sqrt(J);

    successiveDifferenceMean(:, jRank) = differenceMean;
    successiveDifferenceSE(:, jRank) = differenceSE;

    [maximumValue, maximumIndex] = ...
        max(abs(differenceMean));

    maximumSuccessiveDifference(jRank) = maximumValue;
    maximumSuccessiveDifferenceIndex(jRank) = maximumIndex;
    maximumSuccessiveDifferenceSE(jRank) = ...
        differenceSE(maximumIndex);
end

% -------------------------------------------------------------------------
% Paired differences relative to the largest tested rank
% -------------------------------------------------------------------------

referenceRankIndex = numRanks;

referenceDifferenceMean = zeros(Nx, numRanks);
referenceDifferenceSE = zeros(Nx, numRanks);
maximumReferenceDifference = zeros(numRanks, 1);
maximumReferenceDifferenceSE = zeros(numRanks, 1);
maximumReferenceDifferenceIndex = zeros(numRanks, 1);

for jRank = 1:numRanks
    differenceReplicate = ...
        squeeze(FReplicate(:, jRank, :) ...
        - FReplicate(:, referenceRankIndex, :));

    differenceMean = mean(differenceReplicate, 2);
    differenceSE = std(differenceReplicate, 0, 2) / sqrt(J);

    referenceDifferenceMean(:, jRank) = differenceMean;
    referenceDifferenceSE(:, jRank) = differenceSE;

    [maximumValue, maximumIndex] = ...
        max(abs(differenceMean));

    maximumReferenceDifference(jRank) = maximumValue;
    maximumReferenceDifferenceIndex(jRank) = maximumIndex;
    maximumReferenceDifferenceSE(jRank) = ...
        differenceSE(maximumIndex);
end

% -------------------------------------------------------------------------
% Extract matrix diagnostics
% -------------------------------------------------------------------------

relativeFullResidual = zeros(numRanks, 1);
relativeOmittedResidual = zeros(numRanks, 1);
relativeCrossResidual = zeros(numRanks, 1);

for jRank = 1:numRanks
    diagnostics = reductionDiagnostics{jRank};

    relativeFullResidual(jRank) = ...
        diagnostics.relativeMaximumFullResidual;

    relativeOmittedResidual(jRank) = ...
        diagnostics.relativeMaximumOmittedResidual;

    relativeCrossResidual(jRank) = ...
        diagnostics.relativeMaximumCrossResidual;
end

study = struct();

study.rankValues = rankValues;
study.FReplicate = FReplicate;
study.FMean = FMean;
study.FStandardError = FStandardError;

study.pLocalReplicate = pLocalReplicate;
study.pLocalMean = pLocalMean;
study.pLocalStandardError = pLocalStandardError;

study.successiveDifferenceMean = successiveDifferenceMean;
study.successiveDifferenceSE = successiveDifferenceSE;
study.maximumSuccessiveDifference = maximumSuccessiveDifference;
study.maximumSuccessiveDifferenceSE = ...
    maximumSuccessiveDifferenceSE;
study.maximumSuccessiveDifferenceIndex = ...
    maximumSuccessiveDifferenceIndex;

study.referenceDifferenceMean = referenceDifferenceMean;
study.referenceDifferenceSE = referenceDifferenceSE;
study.maximumReferenceDifference = maximumReferenceDifference;
study.maximumReferenceDifferenceSE = ...
    maximumReferenceDifferenceSE;
study.maximumReferenceDifferenceIndex = ...
    maximumReferenceDifferenceIndex;

study.relativeFullResidual = relativeFullResidual;
study.relativeOmittedResidual = relativeOmittedResidual;
study.relativeCrossResidual = relativeCrossResidual;

study.reductionDiagnostics = reductionDiagnostics;
study.referenceRank = rankValues(end);
study.numDirections = K;
study.numReplicates = J;

end

function projected = projectMatrixFamily(matrixFamily, V)

Nx = numel(matrixFamily);
projected = cell(Nx, 1);

for n = 1:Nx
    projected{n} = ...
        V.' * matrixFamily{n} * V;

    projected{n} = ...
        0.5 * (projected{n} + projected{n}.');
end

end
