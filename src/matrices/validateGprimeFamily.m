function diagnostics = validateGprimeFamily(x, G, Gprime)
%VALIDATEGPRIMEFAMILY Validate real-symmetric or complex-Hermitian G'.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx || numel(Gprime) ~= Nx
    error('validateGprimeFamily:LengthMismatch', ...
        'x, G, and Gprime must have matching lengths.');
end

r = size(G{1}, 1);

hermitianDefect = zeros(Nx, 1);
finiteDifferenceDefect = nan(Nx, 1);

for n = 1:Nx
    Gpn = Gprime{n};

    if any(size(Gpn) ~= [r, r])
        error('validateGprimeFamily:DimensionMismatch', ...
            'Gprime{%d} has inconsistent dimensions.', n);
    end

    if any(~isfinite(Gpn(:)))
        error('validateGprimeFamily:NonfiniteValues', ...
            'Gprime{%d} contains nonfinite values.', n);
    end

    hermitianDefect(n) = ...
        norm(Gpn - Gpn', 'fro') ...
        / max(norm(Gpn, 'fro'), eps);
end

for n = 2:(Nx - 1)
    Gfd = ...
        (G{n + 1} - G{n - 1}) ...
        / (x(n + 1) - x(n - 1));

    Gfd = 0.5 * (Gfd + Gfd');

    finiteDifferenceDefect(n) = ...
        norm(Gprime{n} - Gfd, 'fro') ...
        / max(norm(Gprime{n}, 'fro'), eps);
end

diagnostics = struct();
diagnostics.hermitianDefect = hermitianDefect;
diagnostics.finiteDifferenceDefect = finiteDifferenceDefect;
diagnostics.maxHermitianDefect = max(hermitianDefect);

if Nx > 2
    diagnostics.maxFiniteDifferenceDefect = ...
        max(finiteDifferenceDefect(2:end - 1));
else
    diagnostics.maxFiniteDifferenceDefect = NaN;
end

if diagnostics.maxHermitianDefect > 1e-10
    warning('validateGprimeFamily:HermitianDefect', ...
        'Maximum Gprime Hermitian defect is %.3e.', ...
        diagnostics.maxHermitianDefect);
end

% Backward-compatible aliases.
diagnostics.symmetryDefect = hermitianDefect;
diagnostics.maxSymmetryDefect = ...
    diagnostics.maxHermitianDefect;

end
