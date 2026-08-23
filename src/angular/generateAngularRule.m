function [U, weights, meta] = generateAngularRule(r, angularCfg)
%GENERATEANGULARRULE Generate unit-sphere integration directions.
%
% Supported methods:
%
%   'random'
%       Independent Gaussian-normalized directions.
%
%   'circle'
%       Equally spaced projective-circle rule for r = 2.
%
%   'rqmc'
%       Scrambled Sobol' points mapped to the sphere through the
%       inverse-normal transformation and normalization.
%
%   'user'
%       User-supplied directions and weights.
%
% Output:
%   U        r-by-K matrix; each column is a unit direction
%   weights  K-by-1 normalized nonnegative weights
%   meta     angular-rule diagnostics
%
% MATLAB version: R2020b

validateattributes(r, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'r');

if ~isfield(angularCfg, 'method')
    error('generateAngularRule:MissingMethod', ...
        'angularCfg.method is required.');
end

method = lower(char(angularCfg.method));

switch method

    case 'random'
        [U, weights, meta] = generateRandomRule(r, angularCfg);

    case 'circle'
        [U, weights, meta] = generateCircleRule(r, angularCfg);

    case 'rqmc'
        [U, weights, meta] = generateRQMConeReplicate(r, angularCfg);

    case 'user'
        [U, weights, meta] = generateUserRule(r, angularCfg);

    otherwise
        error('generateAngularRule:UnsupportedMethod', ...
            'Unsupported angular method "%s".', ...
            angularCfg.method);
end

% Normalize weights defensively.
weights = weights(:);
weights = weights / sum(weights);

% Normalize all directions defensively.
norms = sqrt(sum(abs(U).^2, 1));

if any(~isfinite(norms)) || any(norms <= 0)
    error('generateAngularRule:InvalidDirectionNorm', ...
        'The angular rule contains a zero or nonfinite direction.');
end

U = U ./ norms;

meta.method = method;
meta.dimension = r;
meta.numDirections = size(U, 2);

end

% =========================================================================
function [U, weights, meta] = generateRandomRule(r, angularCfg)

K = angularCfg.numDirections;

validateattributes(K, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'angularCfg.numDirections');

seed = getOption(angularCfg, 'randomSeed', 1);

rng(seed, 'twister');

if r == 1
    % A quadratic form is invariant under u -> -u, so one projective
    % direction is sufficient. Repeating it K times is harmless.
    U = ones(1, K);
else
    U = randn(r, K);
    U = U ./ sqrt(sum(U.^2, 1));
end

weights = ones(K, 1) / K;

meta = struct();
meta.randomSeed = seed;
meta.isRandomized = true;
meta.replicateIndex = 1;

end

% =========================================================================
function [U, weights, meta] = generateCircleRule(r, angularCfg)

if r ~= 2
    error('generateAngularRule:CircleRequiresR2', ...
        'The circle rule requires stochastic dimension r = 2.');
end

K = angularCfg.numDirections;

validateattributes(K, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'angularCfg.numDirections');

% Integrate over [0,pi) because u and -u give identical quadratic forms.
theta = pi * (0:(K - 1)) / K;

U = [cos(theta); sin(theta)];
weights = ones(K, 1) / K;

meta = struct();
meta.theta = theta(:);
meta.isRandomized = false;
meta.replicateIndex = 1;

end

% =========================================================================
function [U, weights, meta] = generateRQMConeReplicate(r, angularCfg)

if exist('sobolset', 'file') ~= 2
    error('generateAngularRule:SobolUnavailable', ...
        ['sobolset is unavailable. The RQMC implementation requires ', ...
         'the Statistics and Machine Learning Toolbox in MATLAB R2020b.']);
end

K = angularCfg.numDirections;

validateattributes(K, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'angularCfg.numDirections');

seed = getOption(angularCfg, 'randomSeed', 1);
replicateIndex = getOption(angularCfg, 'replicateIndex', 1);

skip = getOption(angularCfg, 'sobolSkip', 1024);
leap = getOption(angularCfg, 'sobolLeap', 0);
applyScramble = getOption(angularCfg, 'scramble', true);
probabilityClip = getOption( ...
    angularCfg, 'probabilityClip', 1e-12);

validateattributes(skip, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'angularCfg.sobolSkip');

validateattributes(leap, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'angularCfg.sobolLeap');

validateattributes(probabilityClip, {'numeric'}, ...
    {'scalar', 'real', 'positive', '<', 0.5}, ...
    mfilename, 'angularCfg.probabilityClip');

% scramble() uses MATLAB's global random stream. Assign a deterministic,
% replicate-specific seed so independent calls are reproducible.
scrambleSeed = seed + 104729 * (replicateIndex - 1);
rng(scrambleSeed, 'twister');

pointSet = sobolset(r, ...
    'Skip', skip, ...
    'Leap', leap);

% Improve lower-dimensional projections where supported.
try
    pointSet = setDigitalShift(pointSet);
catch
    % setDigitalShift is not required and may not be available in all
    % installations. Scrambling below remains the principal randomization.
end

if applyScramble
    pointSet = scramble(pointSet, 'MatousekAffineOwen');
end

% Generate K points in [0,1)^r.
P = net(pointSet, K);

% Protect inverse-normal evaluation from exact 0 and 1.
P = min(max(P, probabilityClip), 1 - probabilityClip);

% Standard-normal inverse:
%
%   Phi^{-1}(p) = -sqrt(2)*erfcinv(2p).
Z = -sqrt(2) * erfcinv(2 * P);

% The rows of Z are sphere-generating vectors. Transpose so columns are
% angular directions.
U = Z.';

norms = sqrt(sum(U.^2, 1));

if any(~isfinite(norms)) || any(norms <= 0)
    error('generateAngularRule:InvalidRQMCDirection', ...
        'RQMC mapping produced a zero or nonfinite direction.');
end

U = U ./ norms;

weights = ones(K, 1) / K;

meta = struct();
meta.isRandomized = true;
meta.replicateIndex = replicateIndex;
meta.scrambleSeed = scrambleSeed;
meta.sobolSkip = skip;
meta.sobolLeap = leap;
meta.scrambled = applyScramble;
meta.probabilityClip = probabilityClip;

end

% =========================================================================
function [U, weights, meta] = generateUserRule(r, angularCfg)

if ~isfield(angularCfg, 'U') || ...
        ~isfield(angularCfg, 'weights')
    error('generateAngularRule:MissingUserRule', ...
        'angularCfg.U and angularCfg.weights are required.');
end

U = angularCfg.U;
weights = angularCfg.weights(:);

if size(U, 1) ~= r
    error('generateAngularRule:UserDimensionMismatch', ...
        'User directions must have r rows.');
end

if numel(weights) ~= size(U, 2)
    error('generateAngularRule:UserWeightMismatch', ...
        'The number of user weights must equal size(U,2).');
end

meta = struct();
meta.isRandomized = false;
meta.replicateIndex = 1;

end

% =========================================================================
function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end

% =========================================================================
function pointSet = setDigitalShift(pointSet)
%SETDIGITALSHIFT Placeholder for compatibility.
%
% No additional transformation is required here. The function exists so
% the call can be safely removed or extended without changing the primary
% RQMC construction.

% Return the point set unchanged.
end