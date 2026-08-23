function Ared = applyComplexReductionToMatrixFamily(A, V)
%APPLYCOMPLEXREDUCTIONTOMATRIXFAMILY Project a Hermitian matrix family.
%
%   Ared{n} = V' * A{n} * V.
%
% MATLAB version: R2020b

Nx = numel(A);
r = size(A{1}, 1);

validateComplexReductionBasis(V, r);

s = size(V, 2);
Ared = cell(Nx, 1);

for n = 1:Nx
    An = A{n};

    if any(size(An) ~= [r, r])
        error('applyComplexReductionToMatrixFamily:DimensionMismatch', ...
            'A{%d} has inconsistent dimensions.', n);
    end

    AredN = V' * An * V;
    AredN = 0.5 * (AredN + AredN');

    if any(size(AredN) ~= [s, s])
        error('applyComplexReductionToMatrixFamily:InternalError', ...
            'Projected matrix has unexpected dimensions.');
    end

    Ared{n} = AredN;
end

end