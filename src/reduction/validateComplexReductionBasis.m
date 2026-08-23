function diagnostics = validateComplexReductionBasis(V, fullRank)
%VALIDATECOMPLEXREDUCTIONBASIS Validate a complex unitary-column basis.
%
% The required condition is:
%
%   V' * V = I.
%
% MATLAB version: R2020b

if ~isnumeric(V) || ndims(V) ~= 2 || isempty(V)
    error('validateComplexReductionBasis:InvalidBasis', ...
        'V must be a nonempty numeric matrix.');
end

if size(V, 1) ~= fullRank
    error('validateComplexReductionBasis:RowDimensionMismatch', ...
        'V has %d rows; expected %d.', ...
        size(V, 1), fullRank);
end

s = size(V, 2);

if s > fullRank
    error('validateComplexReductionBasis:TooManyColumns', ...
        'The reduced rank cannot exceed the full rank.');
end

if any(~isfinite(V(:)))
    error('validateComplexReductionBasis:NonfiniteBasis', ...
        'V contains nonfinite values.');
end

orthogonalityDefect = norm( ...
    V' * V - eye(s), 'fro');

if orthogonalityDefect > 1e-10
    error('validateComplexReductionBasis:OrthogonalityDefect', ...
        'V has unitary-column defect %.6e.', ...
        orthogonalityDefect);
end

diagnostics = struct();
diagnostics.fullRank = fullRank;
diagnostics.reducedRank = s;
diagnostics.orthogonalityDefect = orthogonalityDefect;
diagnostics.isComplex = ~isreal(V);

end
