function study = computeCoupledComplexRankStudy( ...
    x, Gfull, GprimeFull, Vnested, rankValues, ...
    angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTECOUPLEDCOMPLEXRANKSTUDY Coupled proper-complex rank comparison.
%
% One scrambled Sobol' set is generated in 2*rFull real dimensions.
% For rank s, the first s real and first s imaginary normal coordinates
% are combined into complex directions.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);
rFull = size(Gfull{1}, 1);

validateComplexReductionBasis(Vnested, rFull);

rankValues = unique(round(rankValues(:).'));

if any(rankValues < 1) || any(rankValues > rFull)
    error('computeCoupledComplexRankStudy:InvalidRanks', ...
        'All ranks must lie between 1 and rFull.');
end

K = angularCfg.numDirections;
J = angularCfg.numReplicates;
numRanks = numel(rankValues);

weights = ones(K, 1) / K;

FReplicate = zeros(Nx, numRanks, J);
pLocalReplicate = zeros(Nx, numRanks, J);

Gred = cell(numRanks, 1);
GprimeRed = cell(numRanks, 1);

for jRank = 1:numRanks
    s = rankValues(jRank);
    Vs = Vnested(:, 1:s);

    Gred{jRank} = ...
        applyComplexReductionToMatrixFamily(Gfull, Vs);

    if ~isempty(GprimeFull)
        GprimeRed{jRank} = ...
            applyComplexReductionToMatrixFamily( ...
                GprimeFull, Vs);
    else
        GprimeRed{jRank} = [];
    end
end

cfgLocal = struct();
cfgLocal.numerics = numericsCfg;
cfgLocal.maxDetection = maxDetectionCfg;
cfgLocal.gaussianConvention = 'properComplex';

fprintf('\nCoupled proper-complex RQMC rank calculation:\n');
fprintf('  full complex dimension = %d\n', rFull);
fprintf('  directions/replicate   = %d\n', K);
fprintf('  replicates             = %d\n\n', J);

for jRep = 1:J
    fprintf('  replicate %d of %d\n', jRep, J);

    Zfull = generateScrambledSobolNormals( ...
        2 * rFull, K, jRep, angularCfg);

    ZrealFull = Zfull(:, 1:rFull);
    ZimagFull = Zfull(:, (rFull + 1):(2 * rFull));

    for jRank = 1:numRanks
        s = rankValues(jRank);

        Eta = ...
            ZrealFull(:, 1:s) ...
            + 1i * ZimagFull(:, 1:s);

        norms = sqrt(sum(abs(Eta).^2, 2));
        U = bsxfun(@rdivide, Eta.', norms.');

        switch lower(maxDetectionCfg.method)

            case 'grid'
                maxData = computeRunningMaxGrid( ...
                    Gred{jRank}, U, cfgLocal);

            case 'hermite'
                maxData = computeRunningMaxHermite( ...
                    x, ...
                    Gred{jRank}, ...
                    GprimeRed{jRank}, ...
                    U, ...
                    cfgLocal);

            otherwise
                error('computeCoupledComplexRankStudy:UnsupportedMaximum', ...
                    ['The current coupled complex study supports ', ...
                     '''grid'' and ''hermite''.']);
        end

        prob = computeProperComplexTransitionCDF( ...
            maxData.m, weights, s);

        pLocal = ...
            computeProperComplexLocalExceedanceFromGains( ...
                maxData.a, weights, s);

        FReplicate(:, jRank, jRep) = prob.F;
        pLocalReplicate(:, jRank, jRep) = pLocal;
    end
end

FMean = mean(FReplicate, 3);
FStandardError = std(FReplicate, 0, 3) / sqrt(J);

pLocalMean = mean(pLocalReplicate, 3);
pLocalStandardError = ...
    std(pLocalReplicate, 0, 3) / sqrt(J);

referenceRankIndex = numRanks;

referenceDifferenceMean = zeros(Nx, numRanks);
referenceDifferenceSE = zeros(Nx, numRanks);
maximumReferenceDifference = zeros(numRanks, 1);
maximumReferenceDifferenceSE = zeros(numRanks, 1);

for jRank = 1:numRanks
    differenceReplicate = squeeze( ...
        FReplicate(:, jRank, :) ...
        - FReplicate(:, referenceRankIndex, :));

    differenceMean = mean(differenceReplicate, 2);
    differenceSE = ...
        std(differenceReplicate, 0, 2) / sqrt(J);

    referenceDifferenceMean(:, jRank) = differenceMean;
    referenceDifferenceSE(:, jRank) = differenceSE;

    [maximumValue, maximumIndex] = ...
        max(abs(differenceMean));

    maximumReferenceDifference(jRank) = maximumValue;
    maximumReferenceDifferenceSE(jRank) = ...
        differenceSE(maximumIndex);
end

study = struct();

study.rankValues = rankValues;
study.FReplicate = FReplicate;
study.FMean = FMean;
study.FStandardError = FStandardError;

study.pLocalReplicate = pLocalReplicate;
study.pLocalMean = pLocalMean;
study.pLocalStandardError = pLocalStandardError;

study.referenceDifferenceMean = referenceDifferenceMean;
study.referenceDifferenceSE = referenceDifferenceSE;
study.maximumReferenceDifference = ...
    maximumReferenceDifference;
study.maximumReferenceDifferenceSE = ...
    maximumReferenceDifferenceSE;

study.referenceRank = rankValues(end);
study.numDirections = K;
study.numReplicates = J;

end