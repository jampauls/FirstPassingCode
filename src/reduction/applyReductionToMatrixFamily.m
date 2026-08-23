function Ared = applyReductionToMatrixFamily(A, Vred)
%APPLYREDUCTIONTOMATRIXFAMILY Apply a fixed stochastic projection.
%
%   Ared{n} = Vred' * A{n} * Vred

Nx = numel(A);

r = size(A{1}, 1);

if size(Vred, 1) ~= r
    error('applyReductionToMatrixFamily:DimensionMismatch', ...
        'Vred must have r rows.');
end

orthogonalityDefect = norm( ...
    Vred.' * Vred - eye(size(Vred, 2)), ...
    'fro');

if orthogonalityDefect > 1e-10
    warning('applyReductionToMatrixFamily:OrthogonalityDefect', ...
        'Vred has orthogonality defect %.6e.', ...
        orthogonalityDefect);
end

Ared = cell(Nx, 1);

for n = 1:Nx
    Ared{n} = ...
        Vred.' * A{n} * Vred;

    Ared{n} = ...
        0.5 * (Ared{n} + Ared{n}.');
end

end
