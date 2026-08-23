function U = nestedNormalsToDirections(Zfull, rank)
%NESTEDNORMALSTODIRECTIONS Construct an angular rule at nested rank.
%
% Inputs:
%   Zfull  K-by-rFull standard-normal RQMC coordinate array
%   rank   retained stochastic dimension
%
% Output:
%   U      rank-by-K unit directions
%
% MATLAB version: R2020b

if ~isnumeric(Zfull) || ndims(Zfull) ~= 2 || isempty(Zfull)
    error('nestedNormalsToDirections:InvalidInput', ...
        'Zfull must be a nonempty numeric matrix.');
end

validateattributes(rank, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', size(Zfull, 2)}, ...
    mfilename, 'rank');

Z = Zfull(:, 1:rank);

norms = sqrt(sum(Z.^2, 2));

if any(~isfinite(norms)) || any(norms <= 0)
    error('nestedNormalsToDirections:InvalidNorm', ...
        'A nested RQMC direction has zero or nonfinite norm.');
end

U = bsxfun(@rdivide, Z.', norms.');

end
