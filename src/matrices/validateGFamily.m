function diagnostics = validateGFamily(G, cfg)
%VALIDATEGFAMILY Validate a real-symmetric or complex-Hermitian family.
%
% MATLAB version: R2020b

Nx = numel(G);

if Nx == 0
    error('validateGFamily:EmptyInput', ...
        'G must contain at least one matrix.');
end

r = size(G{1}, 1);

diagnostics = struct();
diagnostics.hermitianDefect = zeros(Nx, 1);
diagnostics.minEigenvalue = zeros(Nx, 1);
diagnostics.maxEigenvalue = zeros(Nx, 1);
diagnostics.psdDefect = zeros(Nx, 1);

for n = 1:Nx
    Gn = G{n};

    if ~isnumeric(Gn) || any(size(Gn) ~= [r, r])
        error('validateGFamily:InconsistentSize', ...
            'G{%d} has inconsistent dimensions.', n);
    end

    if any(~isfinite(Gn(:)))
        error('validateGFamily:NonfiniteG', ...
            'G{%d} contains nonfinite values.', n);
    end

    scaleF = max(norm(Gn, 'fro'), eps);

    diagnostics.hermitianDefect(n) = ...
        norm(Gn - Gn', 'fro') / scaleF;

    Gherm = 0.5 * (Gn + Gn');
    lambda = real(eig(Gherm));

    diagnostics.minEigenvalue(n) = min(lambda);
    diagnostics.maxEigenvalue(n) = max(lambda);

    eigScale = max(max(abs(lambda)), eps);

    diagnostics.psdDefect(n) = ...
        max(-diagnostics.minEigenvalue(n), 0) / eigScale;
end

maxHermitianDefect = max(diagnostics.hermitianDefect);
maxPsdDefect = max(diagnostics.psdDefect);

if maxHermitianDefect > 1e-10
    warning('validateGFamily:HermitianDefect', ...
        'Maximum relative Hermitian defect is %.3e.', ...
        maxHermitianDefect);
end

psdTol = 1e-10;

if nargin >= 2 && isfield(cfg, 'numerics') && ...
        isfield(cfg.numerics, 'psdTol')
    psdTol = cfg.numerics.psdTol;
end

if maxPsdDefect > psdTol
    warning('validateGFamily:PsdDefect', ...
        'Maximum relative PSD defect is %.3e.', maxPsdDefect);
end

diagnostics.maxHermitianDefect = maxHermitianDefect;
diagnostics.maxPsdDefect = maxPsdDefect;
diagnostics.isComplex = any(cellfun(@(M) ~isreal(M), G));

% Backward-compatible field name.
diagnostics.symmetryDefect = diagnostics.hermitianDefect;
diagnostics.maxSymmetryDefect = maxHermitianDefect;

end
