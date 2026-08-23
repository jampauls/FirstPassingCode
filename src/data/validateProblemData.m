function validateProblemData(x, B, H, eThresh, cfg)
%VALIDATEPROBLEMDATA Check dimensions and basic problem assumptions.
%
% The factors are assumed to have already been converted to synthesis
% orientation:
%
%   q(x_n) = B{n}*w
%
% where B{n} has dimensions nq-by-r.

if nargin < 5
    cfg = struct();
end

validateattributes(x, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'}, ...
    mfilename, 'x');

x = x(:);
Nx = numel(x);

if any(diff(x) <= 0)
    error('validateProblemData:NonIncreasingX', ...
        'The streamwise stations must be strictly increasing.');
end

if ~iscell(B)
    error('validateProblemData:BNotCell', ...
        'B must be a cell array with one factor per station.');
end

if numel(B) ~= Nx
    error('validateProblemData:BLengthMismatch', ...
        'numel(B) must equal numel(x).');
end

if ~iscell(H)
    error('validateProblemData:HNotCell', ...
        'H must be a cell array with one energy matrix per station.');
end

if numel(H) ~= Nx
    error('validateProblemData:HLengthMismatch', ...
        'numel(H) must equal numel(x).');
end

validateattributes(eThresh, {'numeric'}, ...
    {'vector', 'real', 'finite', 'positive'}, ...
    mfilename, 'eThresh');

eThresh = eThresh(:);

if numel(eThresh) ~= Nx
    error('validateProblemData:ThresholdLengthMismatch', ...
        'numel(eThresh) must equal numel(x).');
end

[nq, r] = size(B{1});

if nq < 1 || r < 1
    error('validateProblemData:EmptyFactor', ...
        'The propagated factor must have nonzero dimensions.');
end

for n = 1:Nx
    Bn = B{n};
    Hn = H{n};

    if ~isnumeric(Bn) || ndims(Bn) ~= 2
        error('validateProblemData:InvalidB', ...
            'B{%d} must be a numeric matrix.', n);
    end

    if any(size(Bn) ~= [nq, r])
        error('validateProblemData:InconsistentBSize', ...
            'B{%d} has size %d-by-%d; expected %d-by-%d.', ...
            n, size(Bn, 1), size(Bn, 2), nq, r);
    end

    if any(~isfinite(Bn(:)))
        error('validateProblemData:NonfiniteB', ...
            'B{%d} contains nonfinite values.', n);
    end

    if ~isnumeric(Hn) || ndims(Hn) ~= 2
        error('validateProblemData:InvalidH', ...
            'H{%d} must be a numeric matrix.', n);
    end

    if any(size(Hn) ~= [nq, nq])
        error('validateProblemData:InconsistentHSize', ...
            'H{%d} has size %d-by-%d; expected %d-by-%d.', ...
            n, size(Hn, 1), size(Hn, 2), nq, nq);
    end

    if any(~isfinite(Hn(:)))
        error('validateProblemData:NonfiniteH', ...
            'H{%d} contains nonfinite values.', n);
    end

    hermitianDefect = norm(Hn - Hn', 'fro') / max(norm(Hn, 'fro'), eps);

    if hermitianDefect > 1e-10
        error('validateProblemData:NonHermitianH', ...
            'H{%d} is not Hermitian. Relative defect = %.3e.', ...
            n, hermitianDefect);
    end
end

if isfield(cfg, 'gaussianConvention')
    if ~strcmpi(cfg.gaussianConvention, 'real')
        warning('validateProblemData:ConventionNotImplemented', ...
            ['The initial code milestone implements the real Gaussian ', ...
             'convention only.']);
    end
end

end