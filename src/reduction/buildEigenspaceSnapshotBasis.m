function [Vsnapshot, info] = ...
    buildEigenspaceSnapshotBasis(G, stationIndices, reductionCfg)
%BUILDEIGENSPACESNAPSHOTBASIS Form a basis from local G eigenspaces.
%
% At each selected station, retain the leading p eigenvectors of G(x).
% The combined snapshot matrix is compressed with an SVD.
%
% MATLAB version: R2020b

Nx = numel(G);
r = size(G{1}, 1);

p = 2;

if isfield(reductionCfg, 'numModesPerStation')
    p = reductionCfg.numModesPerStation;
end

p = min(max(round(p), 1), r);

stationIndices = unique(stationIndices(:));

if any(stationIndices < 1) || ...
        any(stationIndices > Nx)
    error('buildEigenspaceSnapshotBasis:InvalidStation', ...
        'A selected station index is outside the G family.');
end

numStations = numel(stationIndices);
snapshotMatrix = zeros(r, p * numStations);
snapshotEigenvalues = zeros(p, numStations);

column = 0;

for j = 1:numStations
    n = stationIndices(j);

    Gn = 0.5 * (G{n} + G{n}.');

    [V, D] = eig(Gn);
    lambda = real(diag(D));

    [lambda, order] = sort(lambda, 'descend');
    V = real(V(:, order));

    retained = V(:, 1:p);

    snapshotMatrix(:, column + (1:p)) = retained;
    snapshotEigenvalues(:, j) = lambda(1:p);

    column = column + p;
end

[U, S, ~] = svd(snapshotMatrix, 'econ');
singularValues = diag(S);

snapshotTolerance = 1e-12;

if isfield(reductionCfg, 'snapshotTolerance')
    snapshotTolerance = reductionCfg.snapshotTolerance;
end

if isempty(singularValues) || singularValues(1) == 0
    numericalRank = 0;
    Vsnapshot = zeros(r, 0);
else
    numericalRank = nnz( ...
        singularValues ...
        > snapshotTolerance * singularValues(1));

    Vsnapshot = U(:, 1:numericalRank);
end

info = struct();
info.stationIndices = stationIndices;
info.numModesPerStation = p;
info.snapshotEigenvalues = snapshotEigenvalues;
info.singularValues = singularValues;
info.numericalRank = numericalRank;

end
