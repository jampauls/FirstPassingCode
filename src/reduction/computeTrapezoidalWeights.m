function weights = computeTrapezoidalWeights(x)
%COMPUTETRAPEZOIDALWEIGHTS Nonuniform composite trapezoidal weights.
%
% For scalar or matrix-valued samples f(x_n):
%
%   integral f(x) dx approximately equals sum_n weights(n)*f(x_n).
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if Nx < 2
    error('computeTrapezoidalWeights:InsufficientPoints', ...
        'At least two streamwise points are required.');
end

if any(~isfinite(x)) || any(diff(x) <= 0)
    error('computeTrapezoidalWeights:InvalidGrid', ...
        'x must be finite and strictly increasing.');
end

weights = zeros(Nx, 1);

weights(1) = ...
    0.5 * (x(2) - x(1));

weights(end) = ...
    0.5 * (x(end) - x(end - 1));

for n = 2:(Nx - 1)
    weights(n) = ...
        0.5 * (x(n + 1) - x(n - 1));
end

end
