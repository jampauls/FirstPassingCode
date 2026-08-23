function maxData = computeRunningMaxHermite(x, G, Gprime, U, cfg)
%COMPUTERUNNINGMAXHERMITE Resolve directional maxima with Hermite cubics.
%
% The gain for direction u_k is:
%
%   a_k(x) = u_k' G(x) u_k.
%
% Its derivative is:
%
%   a_k'(x) = u_k' Gprime(x) u_k.
%
% A cubic Hermite interpolant is constructed on every streamwise
% interval. Its endpoints and all interior stationary points are
% evaluated to obtain the interval maximum.
%
% Output fields:
%   maxData.a             K-by-Nx endpoint gains
%   maxData.aprime        K-by-Nx endpoint gain derivatives
%   maxData.intervalMax   K-by-(Nx-1) interval maxima
%   maxData.intervalX     K-by-(Nx-1) interval maximizer locations
%   maxData.m             K-by-Nx running maxima through x_n
%   maxData.recordX       K-by-Nx location of active running maximum
%   maxData.info          diagnostics
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx || numel(Gprime) ~= Nx
    error('computeRunningMaxHermite:LengthMismatch', ...
        'x, G, and Gprime must contain the same number of stations.');
end

if Nx < 2
    error('computeRunningMaxHermite:InsufficientStations', ...
        'At least two streamwise stations are required.');
end

[r, K] = size(U);

a = zeros(K, Nx);
aprime = zeros(K, Nx);

% -------------------------------------------------------------------------
% Evaluate endpoint gains and derivatives in batches.
% -------------------------------------------------------------------------

for n = 1:Nx
    if any(size(G{n}) ~= [r, r])
        error('computeRunningMaxHermite:GDimensionMismatch', ...
            'G{%d} has incompatible dimensions.', n);
    end

    if any(size(Gprime{n}) ~= [r, r])
        error('computeRunningMaxHermite:GprimeDimensionMismatch', ...
            'Gprime{%d} has incompatible dimensions.', n);
    end

    Z = G{n} * U;
    gain = real(sum(conj(U) .* Z, 1).');

    Zprime = Gprime{n} * U;
    gainPrime = real(sum(conj(U) .* Zprime, 1).');

    if cfg.numerics.clipSmallNegativeGains
        scale = max(max(abs(gain)), 1);
        smallNegative = gain < 0 & ...
            gain >= -cfg.numerics.psdTol * scale;

        gain(smallNegative) = 0;

        if any(gain < 0)
            warning('computeRunningMaxHermite:NegativeEndpointGain', ...
                ['Material negative directional gain detected at ', ...
                 'station %d.'], n);
        end
    end

    a(:, n) = gain;
    aprime(:, n) = gainPrime;
end

% -------------------------------------------------------------------------
% Hermite interval maxima.
% -------------------------------------------------------------------------

intervalMax = zeros(K, Nx - 1);
intervalX = zeros(K, Nx - 1);

hermiteOptions = struct();
hermiteOptions.polynomialTolerance = 1e-13;
hermiteOptions.rootTolerance = 1e-12;
hermiteOptions.clipNegativeMaximum = true;

totalInteriorCandidates = 0;
totalNegativeClipped = 0;

for n = 1:(Nx - 1)
    [intervalMax(:, n), intervalX(:, n), intervalInfo] = ...
        cubicHermiteIntervalMaximum( ...
            x(n), x(n + 1), ...
            a(:, n), a(:, n + 1), ...
            aprime(:, n), aprime(:, n + 1), ...
            hermiteOptions);

    totalInteriorCandidates = ...
        totalInteriorCandidates ...
        + intervalInfo.totalInteriorCandidates;

    totalNegativeClipped = ...
        totalNegativeClipped ...
        + intervalInfo.numNegativeClipped;
end

% -------------------------------------------------------------------------
% Running maximum through each endpoint.
%
% At x(1), only the inlet value is available. At x(n), n >= 2, include
% the complete interval [x(n-1),x(n)].
% -------------------------------------------------------------------------

m = zeros(K, Nx);
recordX = zeros(K, Nx);

m(:, 1) = max(a(:, 1), 0);
recordX(:, 1) = x(1);

numNewRecords = 0;

for n = 2:Nx
    intervalCandidate = intervalMax(:, n - 1);

    newRecord = intervalCandidate > m(:, n - 1);

    m(:, n) = m(:, n - 1);
    recordX(:, n) = recordX(:, n - 1);

    m(newRecord, n) = intervalCandidate(newRecord);
    recordX(newRecord, n) = intervalX(newRecord, n - 1);

    numNewRecords = numNewRecords + nnz(newRecord);
end

maxData = struct();
maxData.a = a;
maxData.aprime = aprime;
maxData.intervalMax = intervalMax;
maxData.intervalX = intervalX;
maxData.m = m;
maxData.recordX = recordX;

maxData.info = struct();
maxData.info.method = 'hermite';
maxData.info.numDirections = K;
maxData.info.numStations = Nx;
maxData.info.numIntervals = Nx - 1;
maxData.info.totalInteriorCandidates = totalInteriorCandidates;
maxData.info.totalNegativeClipped = totalNegativeClipped;
maxData.info.numNewRecords = numNewRecords;
maxData.info.continuousMaximumResolved = true;
maxData.info.interpolationIsCertifiedBound = false;

end