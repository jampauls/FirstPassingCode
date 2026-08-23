function [Vred, info] = buildTransitionBasis(G, x, cfg)
%BUILDTRANSITIONBASIS Construct a fixed transition-oriented basis.
%
% Supported methods:
%   'integrated'
%   'eigUnion'
%   'hybrid'
%   'greedy'
%
% MATLAB version: R2020b

if ~isfield(cfg, 'reduction')
    error('buildTransitionBasis:MissingConfiguration', ...
        'cfg.reduction is required.');
end

reductionCfg = cfg.reduction;
method = lower(char(reductionCfg.method));

r = size(G{1}, 1);
targetRank = reductionCfg.targetRank;

if isempty(targetRank)
    error('buildTransitionBasis:MissingTargetRank', ...
        'cfg.reduction.targetRank must be specified.');
end

targetRank = min(max(round(targetRank), 1), r);

verbose = true;

if isfield(reductionCfg, 'verbose')
    verbose = reductionCfg.verbose;
end

integratedInfo = struct();
stationInfo = struct();
snapshotInfo = struct();
greedyInfo = struct();

switch method

    case 'integrated'
        [Kint, integratedInfo] = ...
            buildIntegratedTransitionMatrix( ...
                G, x, reductionCfg);

        [V, D] = eig(Kint);
        eigenvalues = real(diag(D));

        [eigenvalues, order] = ...
            sort(eigenvalues, 'descend');

        V = real(V(:, order));

        Vred = V(:, 1:targetRank);

        candidateEigenvalues = eigenvalues;

    case 'eigunion'
        [stationIndices, stationInfo] = ...
            selectTransitionStations( ...
                G, x, reductionCfg);

        [Vsnapshot, snapshotInfo] = ...
            buildEigenspaceSnapshotBasis( ...
                G, stationIndices, reductionCfg);

        if size(Vsnapshot, 2) < targetRank
            warning('buildTransitionBasis:InsufficientSnapshotRank', ...
                ['Snapshot basis rank is %d, below target rank %d. ', ...
                 'Greedy enrichment will be used.'], ...
                size(Vsnapshot, 2), targetRank);
        end

        if size(Vsnapshot, 2) >= targetRank
            Vred = Vsnapshot(:, 1:targetRank);
        else
            [Vred, greedyInfo] = ...
                enrichTransitionBasisGreedy( ...
                    G, Vsnapshot, ...
                    targetRank, reductionCfg);
        end

        candidateEigenvalues = [];

    case 'hybrid'
        [Kint, integratedInfo] = ...
            buildIntegratedTransitionMatrix( ...
                G, x, reductionCfg);

        [Vint, Dint] = eig(Kint);
        lambdaInt = real(diag(Dint));

        [lambdaInt, order] = ...
            sort(lambdaInt, 'descend');

        Vint = real(Vint(:, order));

        numIntegratedModes = min( ...
            reductionCfg.numIntegratedModes, r);

        Vint = Vint(:, 1:numIntegratedModes);

        [stationIndices, stationInfo] = ...
            selectTransitionStations( ...
                G, x, reductionCfg);

        [Vsnapshot, snapshotInfo] = ...
            buildEigenspaceSnapshotBasis( ...
                G, stationIndices, reductionCfg);

        candidateMatrix = [Vint, Vsnapshot];

        [Uc, Sc, ~] = svd(candidateMatrix, 'econ');
        candidateSingularValues = diag(Sc);

        if isempty(candidateSingularValues)
            candidateRank = 0;
        else
            candidateRank = nnz( ...
                candidateSingularValues ...
                > reductionCfg.snapshotTolerance ...
                * candidateSingularValues(1));
        end

        Vcandidate = real(Uc(:, 1:candidateRank));

        if candidateRank >= targetRank
            Vred = Vcandidate(:, 1:targetRank);
        else
            [Vred, greedyInfo] = ...
                enrichTransitionBasisGreedy( ...
                    G, Vcandidate, ...
                    targetRank, reductionCfg);
        end

        candidateEigenvalues = lambdaInt;

    case 'greedy'
        [Vred, greedyInfo] = ...
            enrichTransitionBasisGreedy( ...
                G, zeros(r, 0), ...
                targetRank, reductionCfg);

        candidateEigenvalues = [];

    otherwise
        error('buildTransitionBasis:UnknownMethod', ...
            'Unknown reduction method "%s".', ...
            reductionCfg.method);
end

if size(Vred, 2) < targetRank
    warning('buildTransitionBasis:TargetRankNotReached', ...
        'Constructed rank %d instead of requested rank %d.', ...
        size(Vred, 2), targetRank);
end

% Ensure high-quality orthogonality.
[Vred, ~] = qr(Vred, 0);
Vred = real(Vred);

diagnostics = computeReductionDiagnostics(G, Vred);

info = struct();
info.method = method;
info.fullRank = r;
info.targetRank = targetRank;
info.actualRank = size(Vred, 2);
info.integrated = integratedInfo;
info.stationSelection = stationInfo;
info.snapshot = snapshotInfo;
info.greedy = greedyInfo;
info.matrixDiagnostics = diagnostics;
info.candidateEigenvalues = candidateEigenvalues;

if verbose
    fprintf('Transition-oriented basis:\n');
    fprintf('  method                         = %s\n', method);
    fprintf('  full rank                      = %d\n', r);
    fprintf('  reduced rank                   = %d\n', ...
        size(Vred, 2));
    fprintf('  relative full residual         = %.6e\n', ...
        diagnostics.relativeMaximumFullResidual);
    fprintf('  relative omitted residual      = %.6e\n', ...
        diagnostics.relativeMaximumOmittedResidual);
    fprintf('  relative cross residual        = %.6e\n\n', ...
        diagnostics.relativeMaximumCrossResidual);
end

end