function pLocal = computeLocalExceedance(G, U, weights, r, cfg)
%COMPUTELOCALEXCEEDANCE Compute stationwise threshold exceedance.
%
% The local probability is evaluated using the same radial-directional
% angular rule as the first-transition CDF:
%
%   pLocal(x_n)
%     = E_u[ SF_chi2_r(1 / (u'G(x_n)u)) ].
%
% Using the same U and weights ensures that, for endpoint-only maxima,
%
%   pLocal(x_n) <= F_Xtr(x_n)
%
% direction-by-direction, up to floating-point error.

if nargin < 5
    cfg = struct(); %#ok<NASGU>
end

Nx = numel(G);
[rU, K] = size(U);
weights = weights(:);

if rU ~= r
    error('computeLocalExceedance:DimensionMismatch', ...
        'The row dimension of U must equal r.');
end

if numel(weights) ~= K
    error('computeLocalExceedance:WeightMismatch', ...
        'The number of angular weights must equal size(U,2).');
end

pLocal = zeros(Nx, 1);

for n = 1:Nx
    Gn = G{n};

    Z = Gn * U;
    gain = real(sum(U .* Z, 1).');
    gain(gain < 0 & gain > -1e-12) = 0;

    if any(gain < 0)
        warning('computeLocalExceedance:NegativeGain', ...
            'Material negative directional gain at station %d.', n);
        gain = max(gain, 0);
    end

    positive = gain > 0;
    directionalProbability = zeros(K, 1);

    c = 1 ./ gain(positive);

    directionalProbability(positive) = ...
        gammainc(c / 2, r / 2, 'upper');

    pLocal(n) = weights.' * directionalProbability;
end

pLocal = min(max(pLocal, 0), 1);

end