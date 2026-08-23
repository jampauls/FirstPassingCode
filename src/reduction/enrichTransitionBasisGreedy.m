function [V, info] = ...
    enrichTransitionBasisGreedy(G, Vinitial, targetRank, reductionCfg)
%ENRICHTRANSITIONBASISGREEDY Add largest omitted gain directions.

r = size(G{1}, 1);
Nx = numel(G);

if isempty(Vinitial)
    V = zeros(r, 0);
else
    [V, ~] = qr(Vinitial, 0);
    V = real(V);
end

targetRank = min(targetRank, r);

greedyTolerance = 1e-8;

if isfield(reductionCfg, 'greedyTolerance')
    greedyTolerance = reductionCfg.greedyTolerance;
end

maximumFinalRank = targetRank;

if isfield(reductionCfg, 'maxGreedyModes')
    maximumFinalRank = min( ...
        targetRank, ...
        size(V, 2) + reductionCfg.maxGreedyModes);
end

globalMaximumEigenvalue = 0;

for n = 1:Nx
    Gn = 0.5 * (G{n} + G{n}.');

    globalMaximumEigenvalue = max( ...
        globalMaximumEigenvalue, ...
        max(eig(Gn)));
end

referenceScale = max(globalMaximumEigenvalue, eps);

residualHistory = [];
selectedStationIndices = [];
selectedEigenvalues = [];

while size(V, 2) < targetRank && ...
        size(V, 2) < maximumFinalRank

    if isempty(V)
        Pperp = eye(r);
    else
        Pperp = eye(r) - V * V.';
    end

    bestEigenvalue = -inf;
    bestVector = [];
    bestStation = NaN;

    for n = 1:Nx
        residualMatrix = ...
            Pperp * G{n} * Pperp;

        residualMatrix = ...
            0.5 * (residualMatrix + residualMatrix.');

        [vectors, values] = eig(residualMatrix);
        lambda = real(diag(values));

        [lambdaMaximum, indexMaximum] = max(lambda);

        if lambdaMaximum > bestEigenvalue
            bestEigenvalue = lambdaMaximum;
            bestVector = real(vectors(:, indexMaximum));
            bestStation = n;
        end
    end

    residualHistory(end + 1, 1) = bestEigenvalue; %#ok<AGROW>
    selectedStationIndices(end + 1, 1) = bestStation; %#ok<AGROW>
    selectedEigenvalues(end + 1, 1) = bestEigenvalue; %#ok<AGROW>

    if bestEigenvalue <= ...
            greedyTolerance * referenceScale
        break;
    end

    if ~isempty(V)
        bestVector = ...
            bestVector - V * (V.' * bestVector);

        bestVector = ...
            bestVector - V * (V.' * bestVector);
    end

    vectorNorm = norm(bestVector);

    if vectorNorm <= 1e-12
        warning('enrichTransitionBasisGreedy:DependentVector', ...
            'Greedy enrichment produced a dependent vector.');
        break;
    end

    bestVector = bestVector / vectorNorm;

    V = [V, bestVector]; %#ok<AGROW>

    [V, ~] = qr(V, 0);
    V = real(V);
end

info = struct();
info.residualHistory = residualHistory;
info.selectedStationIndices = selectedStationIndices;
info.selectedEigenvalues = selectedEigenvalues;
info.finalRank = size(V, 2);
info.referenceScale = referenceScale;

end
