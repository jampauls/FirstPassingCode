function G = buildGFromAFamily(A, eThresh, cfg)
%BUILDGFROMAFAMILY Normalize real symmetric or complex Hermitian matrices.

Nx = numel(A);
eThresh = eThresh(:);

if numel(eThresh) ~= Nx
    error('buildGFromAFamily:LengthMismatch', ...
        'eThresh must contain one value per matrix.');
end

r = size(A{1}, 1);
G = cell(Nx, 1);

for n = 1:Nx
    An = A{n};

    if any(size(An) ~= [r, r])
        error('buildGFromAFamily:DimensionMismatch', ...
            'A{%d} has inconsistent dimensions.', n);
    end

    Gn = An / eThresh(n);

    if isfield(cfg, 'numerics') && ...
            isfield(cfg.numerics, 'symmetrizeG') && ...
            cfg.numerics.symmetrizeG

        % Conjugate-transpose symmetrization works for both real symmetric
        % and complex Hermitian matrices.
        Gn = 0.5 * (Gn + Gn');
    end

    G{n} = Gn;
end

end