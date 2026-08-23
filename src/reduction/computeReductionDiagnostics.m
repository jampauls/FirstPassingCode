function diagnostics = ...
    computeReductionDiagnostics(G, V)
%COMPUTEREDUCTIONDIAGNOSTICS Evaluate matrix-family projection errors.
%
% For P = V*V':
%
%   full residual:
%       ||G - PGP||_2
%
%   omitted-subspace residual:
%       ||Pperp*G*Pperp||_2
%
%   cross residual:
%       ||P*G*Pperp||_2
%
% MATLAB version: R2020b

Nx = numel(G);
r = size(G{1}, 1);

validateReductionBasis(V, r);

P = V * V.';
Pperp = eye(r) - P;

fullResidual = zeros(Nx, 1);
omittedResidual = zeros(Nx, 1);
crossResidual = zeros(Nx, 1);
matrixNorm = zeros(Nx, 1);

for n = 1:Nx
    Gn = 0.5 * (G{n} + G{n}.');

    matrixNorm(n) = norm(Gn, 2);

    fullResidual(n) = norm( ...
        Gn - P * Gn * P, 2);

    omittedResidual(n) = norm( ...
        Pperp * Gn * Pperp, 2);

    crossResidual(n) = norm( ...
        P * Gn * Pperp, 2);
end

scale = max(max(matrixNorm), eps);

diagnostics = struct();

diagnostics.matrixNorm = matrixNorm;
diagnostics.fullResidual = fullResidual;
diagnostics.omittedResidual = omittedResidual;
diagnostics.crossResidual = crossResidual;

diagnostics.maximumMatrixNorm = max(matrixNorm);
diagnostics.maximumFullResidual = max(fullResidual);
diagnostics.maximumOmittedResidual = max(omittedResidual);
diagnostics.maximumCrossResidual = max(crossResidual);

diagnostics.relativeMaximumFullResidual = ...
    diagnostics.maximumFullResidual / scale;

diagnostics.relativeMaximumOmittedResidual = ...
    diagnostics.maximumOmittedResidual / scale;

diagnostics.relativeMaximumCrossResidual = ...
    diagnostics.maximumCrossResidual / scale;

end