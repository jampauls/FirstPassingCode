function gain = computeDirectionalGainFamily(matrixFamily, U)
%COMPUTEDIRECTIONALGAINFAMILY Evaluate u.'*A_k*u for all directions.
%
% Inputs:
%   matrixFamily  cell{Nx}, each r-by-r real symmetric
%   U             r-by-K real unit directions
%
% Output:
%   gain          K-by-Nx
%
% MATLAB version: R2020b

Nx = numel(matrixFamily);
[r, K] = size(U);

gain = zeros(K, Nx);

for n = 1:Nx
    An = matrixFamily{n};

    if any(size(An) ~= [r, r])
        error('computeDirectionalGainFamily:DimensionMismatch', ...
            'matrixFamily{%d} is incompatible with U.', n);
    end

    Z = An * U;

    gain(:, n) = real(sum( ...
        U .* Z, 1)).';
end

end
