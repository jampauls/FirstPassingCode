function [Z, meta] = generateScrambledSobolNormals( ...
    dimension, numPoints, replicateIndex, angularCfg)
%GENERATESCRAMPLEDSOBOLNORMALS Generate reusable RQMC normal coordinates.
%
%   Z has size numPoints-by-dimension. For a nested rank s, use:
%
%       Zs = Z(:,1:s);
%       U  = Zs.' ./ sqrt(sum(Zs.^2,2)).';
%
% This couples angular rules across stochastic ranks.
%
% MATLAB version: R2020b

validateattributes(dimension, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'dimension');

validateattributes(numPoints, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'numPoints');

validateattributes(replicateIndex, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'replicateIndex');

if exist('sobolset', 'file') ~= 2
    error('generateScrambledSobolNormals:SobolUnavailable', ...
        'sobolset is unavailable.');
end

baseSeed = getOption(angularCfg, 'randomSeed', 1);
skip = getOption(angularCfg, 'sobolSkip', 1024);
leap = getOption(angularCfg, 'sobolLeap', 0);
applyScramble = getOption(angularCfg, 'scramble', true);
probabilityClip = getOption( ...
    angularCfg, 'probabilityClip', 1e-12);

scrambleSeed = ...
    baseSeed + 104729 * (replicateIndex - 1);

rng(scrambleSeed, 'twister');

pointSet = sobolset(dimension, ...
    'Skip', skip, ...
    'Leap', leap);

if applyScramble
    pointSet = scramble( ...
        pointSet, 'MatousekAffineOwen');
end

P = net(pointSet, numPoints);

P = min(max( ...
    P, probabilityClip), ...
    1 - probabilityClip);

% Inverse standard-normal CDF.
Z = -sqrt(2) * erfcinv(2 * P);

if any(~isfinite(Z(:)))
    error('generateScrambledSobolNormals:NonfiniteOutput', ...
        'The inverse-normal transformation produced nonfinite values.');
end

meta = struct();
meta.dimension = dimension;
meta.numPoints = numPoints;
meta.replicateIndex = replicateIndex;
meta.scrambleSeed = scrambleSeed;
meta.skip = skip;
meta.leap = leap;
meta.scrambled = applyScramble;
meta.probabilityClip = probabilityClip;

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end