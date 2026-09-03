function [Ufamily, diagnostics] = ...
    buildPositiveVariationUpperFamily(A, x, growthReport)
%BUILDPOSITIVEVARIATIONUPPERFAMILY Construct Loewner upper matrices.
%
% If growthReport used:
%
%   D_k = (A_{k+1}-A_k)/dx_k,
%
% then:
%
%   U_1 = A_1,
%   U_{k+1} = U_k + dx_k*D_k^+.
%
% If raw increments were used, the dx factor is omitted.
%
% The resulting sequence satisfies numerically:
%
%   A_j <= U_k for all j <= k,
%
% subject to the positive-part construction tolerance.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(A) ~= Nx
    error('buildPositiveVariationUpperFamily:LengthMismatch', ...
        'A and x must contain the same number of stations.');
end

if numel(growthReport.positivePart) ~= Nx - 1
    error('buildPositiveVariationUpperFamily:IncrementMismatch', ...
        'growthReport does not match the matrix family.');
end

r = size(A{1}, 1);
dx = diff(x);

Ufamily = cell(Nx, 1);

Ucurrent = 0.5 * (A{1} + A{1}.');
Ufamily{1} = Ucurrent;

minimumPairwiseDifferenceEigenvalue = inf(Nx, 1);
minimumPairwiseDifferenceEigenvalue(1) = 0;

minimumAdjacentUpperIncrementEigenvalue = zeros(Nx - 1, 1);

for k = 1:(Nx - 1)
    Dplus = growthReport.positivePart{k};

    if strcmpi(growthReport.incrementScaling, 'perUnitX')
        positiveVariation = dx(k) * Dplus;
    else
        positiveVariation = Dplus;
    end

    positiveVariation = ...
        0.5 * (positiveVariation + positiveVariation.');

    minimumAdjacentUpperIncrementEigenvalue(k) = ...
        min(real(eig(positiveVariation)));

    Ucurrent = Ucurrent + positiveVariation;
    Ucurrent = 0.5 * (Ucurrent + Ucurrent.');

    Ufamily{k + 1} = Ucurrent;

    % Check U_{k+1} - A_j for all j <= k+1.
    minimumDifference = inf;

    for j = 1:(k + 1)
        difference = ...
            Ucurrent ...
            - 0.5 * (A{j} + A{j}.');

        difference = ...
            0.5 * (difference + difference.');

        minimumDifference = min( ...
            minimumDifference, ...
            min(real(eig(difference))));
    end

    minimumPairwiseDifferenceEigenvalue(k + 1) = ...
        minimumDifference;
end

upperScale = max(cellfun( ...
    @(M) norm(M, 2), Ufamily));

diagnostics = struct();

diagnostics.minimumPairwiseDifferenceEigenvalue = ...
    minimumPairwiseDifferenceEigenvalue;

diagnostics.minimumAdjacentUpperIncrementEigenvalue = ...
    minimumAdjacentUpperIncrementEigenvalue;

diagnostics.minimumDominanceEigenvalue = ...
    min(minimumPairwiseDifferenceEigenvalue);

diagnostics.minimumMonotonicIncrementEigenvalue = ...
    min(minimumAdjacentUpperIncrementEigenvalue);

diagnostics.relativeDominanceDefect = ...
    max(-diagnostics.minimumDominanceEigenvalue, 0) ...
    / max(upperScale, eps);

diagnostics.relativeMonotonicityDefect = ...
    max(-diagnostics.minimumMonotonicIncrementEigenvalue, 0) ...
    / max(upperScale, eps);

diagnostics.isNumericallyDominating = ...
    diagnostics.relativeDominanceDefect <= 1e-10;

diagnostics.isNumericallyMonotone = ...
    diagnostics.relativeMonotonicityDefect <= 1e-10;

end
