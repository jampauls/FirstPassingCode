function prob = computeProperComplexTransitionCDF( ...
    m, weights, complexDimension)
%COMPUTEPROPERCOMPLEXTRANSITIONCDF Compute proper-complex first passage.
%
% For:
%
%   z ~ CN(0,I_r),
%
% the squared radius satisfies:
%
%   R^2 ~ Gamma(r,1).
%
% Therefore:
%
%   F_Xtr(x|u)
%       = GammaSF(1/m_u(x); shape=r, scale=1).
%
% MATLAB's normalized upper incomplete gamma function gives:
%
%   gammainc(1/m, r, 'upper').
%
% MATLAB version: R2020b

[K, Nx] = size(m);
weights = weights(:);

if numel(weights) ~= K
    error('computeProperComplexTransitionCDF:WeightMismatch', ...
        'The number of weights must equal size(m,1).');
end

validateattributes(complexDimension, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'complexDimension');

F = zeros(Nx, 1);
S = zeros(Nx, 1);

for n = 1:Nx
    mn = m(:, n);
    positive = mn > 0;

    directionalTransition = zeros(K, 1);
    directionalSurvival = ones(K, 1);

    argument = 1 ./ mn(positive);

    directionalTransition(positive) = ...
        gammainc( ...
            argument, ...
            complexDimension, ...
            'upper');

    directionalSurvival(positive) = ...
        gammainc( ...
            argument, ...
            complexDimension, ...
            'lower');

    Fdirect = ...
        weights.' * directionalTransition;

    Sdirect = ...
        weights.' * directionalSurvival;

    if Fdirect <= Sdirect
        F(n) = Fdirect;
        S(n) = 1 - F(n);
    else
        S(n) = Sdirect;
        F(n) = 1 - S(n);
    end
end

F = min(max(F, 0), 1);
S = min(max(S, 0), 1);

pInterval = zeros(Nx, 1);
pInterval(1) = F(1);

if Nx > 1
    pInterval(2:end) = diff(F);
end

smallNegative = ...
    pInterval < 0 & pInterval > -1e-13;

pInterval(smallNegative) = 0;

if any(pInterval < -1e-12)
    error('computeProperComplexTransitionCDF:NonmonotoneCDF', ...
        'The proper-complex CDF is materially nonmonotone.');
end

prob = struct();

prob.F = F;
prob.S = S;
prob.pInterval = pInterval;
prob.pCensored = S(end);

prob.meta = struct();
prob.meta.gaussianConvention = 'properComplex';
prob.meta.complexDimension = complexDimension;
prob.meta.radialDistribution = 'Gamma(r,1)';

end