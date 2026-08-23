function diagnostics = validateReductionBasis(V, fullRank)
%VALIDATEREDUCTIONBASIS Validate a fixed stochastic reduction basis.
%
%   V must be fullRank-by-s with orthonormal columns.
%
% MATLAB version: R2020b

if ~isnumeric(V) || ndims(V) ~= 2 || isempty(V)
    error('validateReductionBasis:InvalidBasis', ...
        'V must be a nonempty numeric matrix.');
end

if size(V, 1) ~= fullRank
    error('validateReductionBasis:RowDimensionMismatch', ...
        'V has %d rows; expected %d.', ...
        size(V, 1), fullRank);
end

s = size(V, 2);

if s > fullRank
    error('validateReductionBasis:TooManyColumns', ...
        'The reduced rank cannot exceed the full rank.');
end

if any(~isfinite(V(:)))
    error('validateReductionBasis:NonfiniteBasis', ...
        'V contains nonfinite values.');
end

if ~isreal(V)
    error('validateReductionBasis:ComplexBasis', ...
        ['The current real-Gaussian implementation requires a real ', ...
         'stochastic-coordinate basis.']);
end

gramMatrix = V.' * V;
orthogonalityDefect = norm( ...
    gramMatrix - eye(s), 'fro');

if orthogonalityDefect > 1e-10
    error('validateReductionBasis:OrthogonalityDefect', ...
        'V has orthogonality defect %.6e.', ...
        orthogonalityDefect);
end

diagnostics = struct();
diagnostics.fullRank = fullRank;
diagnostics.reducedRank = s;
diagnostics.orthogonalityDefect = orthogonalityDefect;

end