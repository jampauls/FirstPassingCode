function maxData = computeRunningMaxGrid(G, U, cfg)
%COMPUTERUNNINGMAXGRID Compute endpoint-only directional running maxima.
%
% Inputs:
%   G   cell array of r-by-r normalized energy matrices
%   U   r-by-K angular direction matrix
%
% Output fields:
%   maxData.a       K-by-Nx endpoint directional gains
%   maxData.m       K-by-Nx running maxima
%   maxData.info    diagnostic structure

Nx = numel(G);
[r, K] = size(U);

a = zeros(K, Nx);
m = zeros(K, Nx);

numNegativeClipped = 0;
minimumRawGain = inf;

for n = 1:Nx
    Gn = G{n};

    if any(size(Gn) ~= [r, r])
        error('computeRunningMaxGrid:DimensionMismatch', ...
            'G{%d} is not compatible with U.', n);
    end

    % Z(:,k) = G*u_k.
    Z = Gn * U;

    % Columnwise quadratic forms u_k' G u_k.
    gain = sum(conj(U) .* Z, 1).';
    gain = real(gain);

    minimumRawGain = min(minimumRawGain, min(gain));

    if isfield(cfg, 'numerics') && ...
            isfield(cfg.numerics, 'clipSmallNegativeGains') && ...
            cfg.numerics.clipSmallNegativeGains

        negative = gain < 0;
        numNegativeClipped = numNegativeClipped + nnz(negative);
        gain(negative) = 0;
    end

    a(:, n) = gain;

    if n == 1
        m(:, n) = gain;
    else
        m(:, n) = max(m(:, n - 1), gain);
    end
end

maxData = struct();
maxData.a = a;
maxData.m = m;

maxData.info = struct();
maxData.info.method = 'grid';
maxData.info.numDirections = K;
maxData.info.numStations = Nx;
maxData.info.numNegativeClipped = numNegativeClipped;
maxData.info.minimumRawGain = minimumRawGain;
maxData.info.continuousMaximumResolved = false;

end