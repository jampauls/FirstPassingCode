function weights = localPolynomialDerivativeWeights( ...
    evaluationPoint, stencilPoints, derivativeOrder)
%LOCALPOLYNOMIALDERIVATIVEWEIGHTS Local polynomial FD weights.
%
%   weights approximate the derivative according to:
%
%       f^(m)(x0) approximately equals sum_j weights(j)*f(x_j).
%
% The weights differentiate polynomials through degree n-1 exactly,
% subject to floating-point accuracy.
%
% MATLAB version: R2020b

if nargin < 3
    derivativeOrder = 1;
end

stencilPoints = stencilPoints(:);
numPoints = numel(stencilPoints);

validateattributes(evaluationPoint, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, ...
    mfilename, 'evaluationPoint');

validateattributes(derivativeOrder, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', '<', numPoints}, ...
    mfilename, 'derivativeOrder');

if any(~isfinite(stencilPoints)) || ...
        numel(unique(stencilPoints)) ~= numPoints
    error('localPolynomialDerivativeWeights:InvalidStencil', ...
        'Stencil points must be finite and distinct.');
end

dx = stencilPoints - evaluationPoint;
scale = max(abs(dx));

if scale == 0
    error('localPolynomialDerivativeWeights:ZeroScale', ...
        'The differentiation stencil has zero scale.');
end

scaledDx = dx / scale;

% Moment equations:
%
%   sum_j c_j * scaledDx_j^p
%       = factorial(p), if p = derivativeOrder
%       = 0,            otherwise.
momentMatrix = zeros(numPoints, numPoints);

for p = 0:(numPoints - 1)
    momentMatrix(p + 1, :) = ...
        scaledDx.'.^p;
end

rightHandSide = zeros(numPoints, 1);
rightHandSide(derivativeOrder + 1) = ...
    factorial(derivativeOrder);

scaledWeights = ...
    momentMatrix \ rightHandSide;

weights = ...
    scaledWeights / scale^derivativeOrder;

end