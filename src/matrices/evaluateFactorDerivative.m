function Bprime = evaluateFactorDerivative( ...
    factorArray, x, stationIndex, stencilSize)
%EVALUATEFACTORDERIVATIVE Differentiate a propagated factor in x.
%
% Inputs:
%   factorArray   physicalDimension-by-Nx-by-r array
%   x             Nx-by-1 streamwise coordinate
%   stationIndex  derivative evaluation station
%   stencilSize   number of local polynomial points
%
% Output:
%   Bprime        physicalDimension-by-r derivative matrix
%
% MATLAB version: R2020b

x = x(:);

[physicalDimension, Nx, r] = size(factorArray);

if numel(x) ~= Nx
    error('evaluateFactorDerivative:CoordinateLengthMismatch', ...
        'x and factorArray have inconsistent station counts.');
end

indices = selectLocalStencil( ...
    Nx, stationIndex, stencilSize);

weights = localPolynomialDerivativeWeights( ...
    x(stationIndex), ...
    x(indices), ...
    1);

Bprime = complex(zeros(physicalDimension, r));

for j = 1:numel(indices)
    Bj = squeeze(factorArray(:, indices(j), :));

    Bprime = Bprime + weights(j) * Bj;
end

end
