function G = buildGFamily(B, H, eThresh, cfg)
%BUILDGFAMILY Construct normalized stochastic-space energy matrices.
%
% For real Gaussian coordinates and possibly complex physical states:
%
%   A(x) = real(B(x)' * H(x) * B(x))
%   G(x) = A(x) / eThresh(x)
%
% Inputs:
%   B       cell array, B{n} is nq-by-r
%   H       cell array, H{n} is nq-by-nq
%   eThresh Nx-by-1 positive threshold values
%   cfg     configuration structure
%
% Output:
%   G       cell array, G{n} is r-by-r real symmetric

Nx = numel(B);
eThresh = eThresh(:);

if numel(H) ~= Nx || numel(eThresh) ~= Nx
    error('buildGFamily:InputLengthMismatch', ...
        'B, H, and eThresh must contain the same number of stations.');
end

G = cell(Nx, 1);

for n = 1:Nx
    Bn = B{n};
    Hn = H{n};

    An = real(Bn' * Hn * Bn);

    if isfield(cfg, 'numerics') && ...
            isfield(cfg.numerics, 'symmetrizeG') && ...
            cfg.numerics.symmetrizeG
        An = 0.5 * (An + An.');
    end

    Gn = An / eThresh(n);

    if isfield(cfg, 'numerics') && ...
            isfield(cfg.numerics, 'psdProjection') && ...
            cfg.numerics.psdProjection

        [V, D] = eig(0.5 * (Gn + Gn.'));
        lambda = real(diag(D));

        scale = max(max(abs(lambda)), 1);
        tol = cfg.numerics.psdTol * scale;

        if min(lambda) < -tol
            warning('buildGFamily:MaterialNegativeEigenvalue', ...
                ['G{%d} has a negative eigenvalue %.3e that exceeds ', ...
                 'the PSD tolerance.'], n, min(lambda));
        end

        lambda(lambda < 0) = 0;
        Gn = V * diag(lambda) * V.';
        Gn = real(0.5 * (Gn + Gn.'));
    end

    G{n} = Gn;
end

end
