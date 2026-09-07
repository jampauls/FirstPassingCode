function [Gred, Vred, info] = reduceEnsembleMode(G, x, perModeCfg)
%REDUCEENSEMBLEMODE Reduce one ensemble mode to a small stochastic rank.
%
%   [Gred, Vred, info] = reduceEnsembleMode(G, x, perModeCfg)
%
% G is the cell array of this mode's normalized energy matrices
% G_i(x) = A_i(x)/eThreshTotal(x). perModeCfg reuses the same fields as
% cfg.reduction (method, targetRank, etc.), applied independently to this
% mode.
%
% MATLAB version: R2020b

r = size(G{1}, 1);

method = lower(char(perModeCfg.method));

if strcmp(method, 'none') || perModeCfg.targetRank >= r
    Gred = G;
    Vred = eye(r);

    info = struct();
    info.method = 'none';
    info.targetRank = r;
    info.fullRank = r;

    return;
end

modeCfg = struct();
modeCfg.reduction = perModeCfg;

[Vred, info] = buildTransitionBasis(G, x, modeCfg);

Gred = applyReductionToMatrixFamily(G, Vred);

info.fullRank = r;
info.targetRank = size(Vred, 2);

end
