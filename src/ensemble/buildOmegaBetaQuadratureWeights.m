function [weights, info] = buildOmegaBetaQuadratureWeights(omega, beta, cfg)
%BUILDOMEGABETAQUADRATUREWEIGHTS Nearest-cell weights on an irregular grid.

omega = omega(:);
beta = beta(:);
M = numel(omega);

if numel(beta) ~= M || M < 3 || any(~isfinite(omega)) || ...
        any(~isfinite(beta))
    error('buildOmegaBetaQuadratureWeights:InvalidPoints', ...
        'omega and beta must contain at least 3 finite paired values.');
end

points = [omega, beta];
if size(unique(points, 'rows'), 1) ~= M
    error('buildOmegaBetaQuadratureWeights:DuplicatePoints', ...
        'Each omega-beta pair must occur exactly once.');
end

omegaLimits = [min(omega), 0.1];
betaLimits = [0, 1];
numOmegaCells = 500;
numBetaCells = 500;

if isfield(cfg, 'ensembleQuadrature')
    quadratureCfg = cfg.ensembleQuadrature;
    if isfield(quadratureCfg, 'omegaLimits') && ...
            ~isempty(quadratureCfg.omegaLimits)
        omegaLimits = quadratureCfg.omegaLimits;
    end
    if isfield(quadratureCfg, 'betaLimits') && ...
            ~isempty(quadratureCfg.betaLimits)
        betaLimits = quadratureCfg.betaLimits;
    end
    if isfield(quadratureCfg, 'numOmegaCells') && ...
            ~isempty(quadratureCfg.numOmegaCells)
        numOmegaCells = quadratureCfg.numOmegaCells;
    end
    if isfield(quadratureCfg, 'numBetaCells') && ...
            ~isempty(quadratureCfg.numBetaCells)
        numBetaCells = quadratureCfg.numBetaCells;
    end
end

validateattributes(omegaLimits, {'numeric'}, {'vector', 'numel', 2, 'finite'});
validateattributes(betaLimits, {'numeric'}, {'vector', 'numel', 2, 'finite'});
validateattributes(numOmegaCells, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(numBetaCells, {'numeric'}, {'scalar', 'integer', '>=', 2});

if omegaLimits(2) <= omegaLimits(1) || betaLimits(2) <= betaLimits(1)
    error('buildOmegaBetaQuadratureWeights:InvalidLimits', ...
        'Quadrature limits must be strictly increasing.');
end

dOmega = diff(omegaLimits) / numOmegaCells;
dBeta = diff(betaLimits) / numBetaCells;
omegaCenters = omegaLimits(1) + dOmega * ((0:(numOmegaCells - 1)) + 0.5);
betaCenters = betaLimits(1) + dBeta * ((0:(numBetaCells - 1)) + 0.5);

weights = zeros(M, 1);
hullWeights = zeros(M, 1);
hullIndex = convhull(omega, beta);

for j = 1:numBetaCells
    pointOmega = omegaCenters(:);
    pointBeta = betaCenters(j);
    distanceSquared = (pointOmega - omega.').^2 + (pointBeta - beta.').^2;
    [~, nearest] = min(distanceSquared, [], 2);
    weights = weights + accumarray(nearest, dOmega * dBeta, [M, 1]);

    insideHull = inpolygon(pointOmega, pointBeta * ones(numOmegaCells, 1), ...
        omega(hullIndex), beta(hullIndex));
    hullWeights = hullWeights + accumarray( ...
        nearest(insideHull), dOmega * dBeta, [M, 1]);
end

domainArea = diff(omegaLimits) * diff(betaLimits);

info = struct();
info.method = 'nearestCell';
info.omegaLimits = omegaLimits(:).';
info.betaLimits = betaLimits(:).';
info.numOmegaCells = numOmegaCells;
info.numBetaCells = numBetaCells;
info.domainArea = domainArea;
info.weightSum = sum(weights);
info.hullWeightSum = sum(hullWeights);
info.boundaryExtrapolationFraction = 1 - sum(hullWeights) / domainArea;
info.weightsZeroOutsideHull = hullWeights;
info.weightL1DifferenceZeroOutsideHull = sum(abs(weights - hullWeights));

end