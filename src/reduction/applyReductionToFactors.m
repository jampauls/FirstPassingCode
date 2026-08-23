function Bred = applyReductionToFactors(B, Vred)
%APPLYREDUCTIONTOFACTORS Apply a fixed stochastic-coordinate basis.

Nx = numel(B);
Bred = cell(Nx, 1);

if size(B{1}, 2) ~= size(Vred, 1)
    error('applyReductionToFactors:DimensionMismatch', ...
        'The row dimension of Vred must equal the stochastic rank.');
end

orthogonalityDefect = ...
    norm(Vred.' * Vred - eye(size(Vred, 2)), 'fro');

if orthogonalityDefect > 1e-10
    warning('applyReductionToFactors:NonorthogonalBasis', ...
        'Vred has orthogonality defect %.3e.', ...
        orthogonalityDefect);
end

for n = 1:Nx
    Bred{n} = B{n} * Vred;
end

end