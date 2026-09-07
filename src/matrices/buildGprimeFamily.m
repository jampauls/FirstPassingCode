function [Gprime, info] = buildGprimeFamily(x, G, meta, cfg, Vred)
%BUILDGPRIMEFAMILY Construct streamwise derivatives G'(x).
%
%   [Gprime, info] = buildGprimeFamily(x, G, meta, cfg, Vred)
%
% Priority:
%
%   1. Use meta.GprimeExact when available.
%   2. Otherwise use local three-point polynomial differentiation of G.
%
% If Vred is supplied and the exact derivatives are in the unreduced
% stochastic coordinates, they are projected as:
%
%   GprimeReduced = Vred' * GprimeFull * Vred.
%
% Inputs:
%   x       Nx-by-1 strictly increasing streamwise grid
%   G       cell array of matrices used by the probability calculation
%   meta    metadata structure
%   cfg     configuration structure
%   Vred    optional fixed reduction basis
%
% Outputs:
%   Gprime  cell array of matrices matching G
%   info    derivative-method diagnostics
%
% MATLAB version: R2020b

if nargin < 5
    Vred = [];
end

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx
    error('buildGprimeFamily:LengthMismatch', ...
        'numel(G) must equal numel(x).');
end

if Nx < 2
    error('buildGprimeFamily:InsufficientStations', ...
        'At least two streamwise stations are required.');
end

r = size(G{1}, 1);

info = struct();
info.usedExactDerivative = false;
info.usedFiniteDifference = false;
info.method = '';

% -------------------------------------------------------------------------
% Use exact or supplied derivative data when available.
% -------------------------------------------------------------------------

if isfield(meta, 'GprimeExact') && ...
        ~isempty(meta.GprimeExact)

    GprimeSource = meta.GprimeExact;

    if numel(GprimeSource) ~= Nx
        error('buildGprimeFamily:ExactLengthMismatch', ...
            'meta.GprimeExact must contain one matrix per station.');
    end

    sourceRank = size(GprimeSource{1}, 1);

    if sourceRank == r
        % Exact data already match the matrices used in the probability
        % calculation.
        Gprime = GprimeSource;

    elseif ~isempty(Vred) && ...
            sourceRank == size(Vred, 1) && ...
            r == size(Vred, 2)

        Gprime = cell(Nx, 1);

        for n = 1:Nx
            Gprime{n} = ...
                Vred.' * GprimeSource{n} * Vred;
        end

    else
        error('buildGprimeFamily:ExactDimensionMismatch', ...
            ['Exact Gprime matrices do not match the active stochastic ', ...
             'dimension and cannot be projected using Vred.']);
    end

    info.usedExactDerivative = true;
    info.method = 'metadataExact';

else
    % ---------------------------------------------------------------------
    % Finite-difference fallback.
    % ---------------------------------------------------------------------

    Gprime = differentiateMatrixFamilyThreePoint(x, G);

    info.usedFiniteDifference = true;
    info.method = 'threePointPolynomial';
end

% -------------------------------------------------------------------------
% Symmetry and dimension cleanup.
% -------------------------------------------------------------------------

for n = 1:Nx
    Gpn = Gprime{n};

    if any(size(Gpn) ~= [r, r])
        error('buildGprimeFamily:OutputDimensionMismatch', ...
            'Gprime{%d} does not match the active G dimension.', n);
    end

    % Preserve complex Hermitian derivatives.
    Gpn = 0.5 * (Gpn + Gpn');

    if any(~isfinite(Gpn(:)))
        error('buildGprimeFamily:NonfiniteDerivative', ...
            'Gprime{%d} contains nonfinite values.', n);
    end

    Gprime{n} = Gpn;
end

info.numStations = Nx;
info.dimension = r;
info.cfg = cfg;

end

function Gprime = differentiateMatrixFamilyThreePoint(x, G)
%DIFFERENTIATEMATRIXFAMILYTHREEPOINT Differentiate matrix data using
% local degree-two polynomial weights on a nonuniform grid.

Nx = numel(x);
r = size(G{1}, 1);

Gprime = cell(Nx, 1);

if Nx == 2
    slope = (G{2} - G{1}) / (x(2) - x(1));
    Gprime{1} = slope;
    Gprime{2} = slope;
    return;
end

for n = 1:Nx
    if n == 1
        indices = 1:3;
    elseif n == Nx
        indices = (Nx - 2):Nx;
    else
        indices = (n - 1):(n + 1);
    end

    xCenter = x(n);
    dx = x(indices) - xCenter;

    % Choose weights c_j such that:
    %
    %   sum_j c_j            = 0
    %   sum_j c_j*dx_j       = 1
    %   sum_j c_j*dx_j^2     = 0
    %
    % This differentiates every quadratic polynomial exactly at xCenter.
    momentMatrix = [ ...
        ones(1, 3); ...
        dx.'; ...
        (dx.^2).' ...
    ];

    weights = momentMatrix \ [0; 1; 0];

    % A sparse accumulator preserves block-diagonal sparsity when G is
    % sparse (e.g. the ensemble workflow), while remaining correct for
    % dense G since sparse-plus-dense addition yields a dense result.
    if issparse(G{indices(1)})
        Gpn = sparse(r, r);
    else
        Gpn = zeros(r, r);
    end

    for j = 1:3
        Gpn = Gpn + weights(j) * G{indices(j)};
    end

    Gprime{n} = 0.5 * (Gpn + Gpn');
end

end
