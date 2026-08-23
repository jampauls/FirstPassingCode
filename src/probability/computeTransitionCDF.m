function prob = computeTransitionCDF(m, weights, r, cfg)
%COMPUTETRANSITIONCDF Compute survival and first-transition probabilities.
%
% For real Gaussian coordinates:
%
%   S(x_n) = sum_k weights(k) * F_chi2_r(1/m_k(x_n))
%
%   F(x_n) = sum_k weights(k) * SF_chi2_r(1/m_k(x_n))
%
% The implementation uses gammainc, so the Statistics and Machine
% Learning Toolbox is not required.
%
% MATLAB identities:
%
%   F_chi2_r(z)  = gammainc(z/2, r/2, 'lower')
%   SF_chi2_r(z) = gammainc(z/2, r/2, 'upper')

if nargin < 4
    cfg = struct();
end

[K, Nx] = size(m);
weights = weights(:);

if numel(weights) ~= K
    error('computeTransitionCDF:WeightLengthMismatch', ...
        'The number of weights must equal size(m,1).');
end

if any(weights < 0) || abs(sum(weights) - 1) > 1e-12
    error('computeTransitionCDF:InvalidWeights', ...
        'Weights must be nonnegative and sum to one.');
end

validateattributes(r, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'r');

S = zeros(Nx, 1);
F = zeros(Nx, 1);

for n = 1:Nx
    mn = m(:, n);

    c = inf(K, 1);
    positive = mn > 0;
    c(positive) = 1 ./ mn(positive);

    directionalSurvival = ones(K, 1);
    directionalTransition = zeros(K, 1);

    % chi-square CDF and survival function.
    directionalSurvival(positive) = ...
        gammainc(c(positive) / 2, r / 2, 'lower');

    directionalTransition(positive) = ...
        gammainc(c(positive) / 2, r / 2, 'upper');

    % Sum both tails independently for diagnostics.
    Sdirect = weights.' * directionalSurvival;
    Fdirect = weights.' * directionalTransition;

    % Choose the smaller tail as the directly accumulated quantity and
    % obtain the larger by complementarity.
    if Fdirect <= Sdirect
        F(n) = Fdirect;
        S(n) = 1 - F(n);
    else
        S(n) = Sdirect;
        F(n) = 1 - S(n);
    end
end

% Remove roundoff-scale excursions.
S = min(max(S, 0), 1);
F = min(max(F, 0), 1);

% Inlet mass and interval transition probabilities.
pInterval = zeros(Nx, 1);
pInterval(1) = F(1);

if Nx > 1
    pInterval(2:end) = diff(F);
end

smallNegative = pInterval < 0 & pInterval > -1e-12;
pInterval(smallNegative) = 0;

pCensored = S(end);

% Discrete conditional transition probability.
hDiscrete = nan(Nx, 1);
hDiscrete(1) = F(1);

for n = 2:Nx
    if S(n - 1) > 0
        hDiscrete(n) = ...
            (S(n - 1) - S(n)) / S(n - 1);
    end
end

prob = struct();
prob.S = S;
prob.F = F;
prob.pInterval = pInterval;
prob.pCensored = pCensored;
prob.hDiscrete = hDiscrete;

prob.meta = struct();
prob.meta.gaussianConvention = 'real';
prob.meta.radialDegreesOfFreedom = r;
prob.meta.numDirections = K;
prob.meta.totalProbability = sum(pInterval) + pCensored;
prob.meta.cfg = cfg;

end