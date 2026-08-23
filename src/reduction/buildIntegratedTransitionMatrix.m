function [Kint, info] = ...
    buildIntegratedTransitionMatrix(G, x, reductionCfg)
%BUILDINTEGRATEDTRANSITIONMATRIX Build a finite-horizon aggregate matrix.
%
%   Kint = integral omega(x)*G(x) dx
%
% Supported normalization:
%
%   'uniform'
%       Integrate G(x) directly.
%
%   'traceNormalized'
%       Integrate G(x)/trace(G(x)).
%
%   'spectralNormalized'
%       Integrate G(x)/lambda_max(G(x)).
%
% The normalized variants emphasize directional subspace coverage rather
% than stations with the largest absolute gain.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx
    error('buildIntegratedTransitionMatrix:LengthMismatch', ...
        'G must contain one matrix per streamwise station.');
end

r = size(G{1}, 1);

integrationWeights = computeTrapezoidalWeights(x);

normalization = 'uniform';

if isfield(reductionCfg, 'integrationWeight') && ...
        ~isempty(reductionCfg.integrationWeight)
    normalization = ...
        lower(char(reductionCfg.integrationWeight));
end

Kint = zeros(r, r);
normalizationScale = ones(Nx, 1);

for n = 1:Nx
    Gn = 0.5 * (G{n} + G{n}.');

    switch normalization

        case 'uniform'
            scale = 1;

        case 'tracenormalized'
            scale = trace(Gn);

        case 'spectralnormalized'
            scale = max(eig(Gn));

        otherwise
            error('buildIntegratedTransitionMatrix:UnknownWeight', ...
                'Unknown integrationWeight "%s".', ...
                reductionCfg.integrationWeight);
    end

    if scale <= eps
        normalizationScale(n) = 0;
        continue;
    end

    normalizationScale(n) = scale;

    Kint = Kint ...
        + integrationWeights(n) * Gn / scale;
end

Kint = 0.5 * (Kint + Kint.');

info = struct();
info.integrationWeights = integrationWeights;
info.normalization = normalization;
info.normalizationScale = normalizationScale;
info.trace = trace(Kint);
info.maximumEigenvalue = max(eig(Kint));

end