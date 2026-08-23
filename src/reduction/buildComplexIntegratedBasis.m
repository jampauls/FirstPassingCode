function [V, info] = buildComplexIntegratedBasis( ...
    G, x, targetRank, options)
%BUILDCOMPLEXINTEGRATEDBASIS Build a complex transition-oriented basis.
%
%   K = integral omega(x) G(x) dx
%
% where G(x) is Hermitian. The leading eigenvectors of K define a fixed
% complex stochastic-coordinate basis.
%
% Supported normalization:
%   'uniform'
%   'traceNormalized'
%   'spectralNormalized'
%
% MATLAB version: R2020b

if nargin < 4
    options = struct();
end

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx
    error('buildComplexIntegratedBasis:LengthMismatch', ...
        'G must contain one matrix per streamwise station.');
end

r = size(G{1}, 1);

validateattributes(targetRank, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', r}, ...
    mfilename, 'targetRank');

normalization = getOption( ...
    options, 'integrationWeight', 'uniform');

verbose = getOption(options, 'verbose', true);

integrationWeights = computeTrapezoidalWeights(x);

Kint = complex(zeros(r, r));
normalizationScale = ones(Nx, 1);

for n = 1:Nx
    Gn = 0.5 * (G{n} + G{n}');

    switch lower(char(normalization))

        case 'uniform'
            scale = 1;

        case 'tracenormalized'
            scale = real(trace(Gn));

        case 'spectralnormalized'
            scale = max(real(eig(Gn)));

        otherwise
            error('buildComplexIntegratedBasis:UnknownWeight', ...
                'Unknown integrationWeight "%s".', normalization);
    end

    if scale <= eps
        normalizationScale(n) = 0;
        continue;
    end

    normalizationScale(n) = scale;

    Kint = Kint ...
        + integrationWeights(n) * Gn / scale;
end

Kint = 0.5 * (Kint + Kint');

[Vall, Dall] = eig(Kint);
lambda = real(diag(Dall));

[lambda, order] = sort(lambda, 'descend');
Vall = Vall(:, order);

V = Vall(:, 1:targetRank);

% Fix arbitrary complex phases for reproducibility.
for j = 1:targetRank
    [~, indexLargest] = max(abs(V(:, j)));
    phase = angle(V(indexLargest, j));
    V(:, j) = V(:, j) * exp(-1i * phase);

    if real(V(indexLargest, j)) < 0
        V(:, j) = -V(:, j);
    end
end

validateComplexReductionBasis(V, r);

positiveLambda = max(lambda, 0);

if sum(positiveLambda) > 0
    cumulativeFraction = ...
        cumsum(positiveLambda) / sum(positiveLambda);
else
    cumulativeFraction = zeros(size(lambda));
end

info = struct();
info.method = 'complexIntegrated';
info.targetRank = targetRank;
info.fullRank = r;
info.integrationWeight = normalization;
info.integrationWeights = integrationWeights;
info.normalizationScale = normalizationScale;
info.integratedMatrix = Kint;
info.eigenvalues = lambda;
info.cumulativeFraction = cumulativeFraction;
info.retainedFraction = cumulativeFraction(targetRank);

if verbose
    fprintf('Complex integrated basis:\n');
    fprintf('  full rank          = %d\n', r);
    fprintf('  reduced rank       = %d\n', targetRank);
    fprintf('  retained fraction  = %.10f\n\n', ...
        info.retainedFraction);
end

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end