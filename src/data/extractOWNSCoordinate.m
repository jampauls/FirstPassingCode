function [xCoordinate, coordinateMeta] = ...
    extractOWNSCoordinate(solution, coordinateCfg)
%EXTRACTOWNSCOORDINATE Extract one ordered first-passage coordinate.
%
% Supported methods:
%
%   'wallArcLength'
%       Three-dimensional arc length along wall-normal row 1.
%
%   'referenceArcLength'
%       Arc length along coordinateCfg.referenceIndex.
%
%   'wallX'
%       Physical x coordinate along wall-normal row 1.
%
%   'referenceX'
%       Physical x coordinate along coordinateCfg.referenceIndex.
%
%   'xi'
%       OWNS marching coordinate solution.xi.
%
% MATLAB version: R2020b

if nargin < 2
    coordinateCfg = struct();
end

if ~isfield(coordinateCfg, 'method')
    coordinateCfg.method = 'wallArcLength';
end

if ~isfield(coordinateCfg, 'referenceIndex')
    coordinateCfg.referenceIndex = 1;
end

if ~isfield(coordinateCfg, 'includeZ')
    coordinateCfg.includeZ = true;
end

if ~isfield(coordinateCfg, 'zeroOrigin')
    coordinateCfg.zeroOrigin = true;
end

method = lower(char(coordinateCfg.method));

Ny = size(solution.x, 1);
Nx = size(solution.x, 2);

referenceIndex = coordinateCfg.referenceIndex;

validateattributes(referenceIndex, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, '<=', Ny}, ...
    mfilename, 'coordinateCfg.referenceIndex');

switch method

    case 'wallarclength'
        referenceIndex = 1;
        [xCoordinate, physicalCoordinate] = ...
            computeArcLength(solution, referenceIndex, ...
                coordinateCfg.includeZ);

    case 'referencearclength'
        [xCoordinate, physicalCoordinate] = ...
            computeArcLength(solution, referenceIndex, ...
                coordinateCfg.includeZ);

    case 'wallx'
        referenceIndex = 1;
        physicalCoordinate = solution.x(referenceIndex, :).';
        xCoordinate = physicalCoordinate;

    case 'referencex'
        physicalCoordinate = solution.x(referenceIndex, :).';
        xCoordinate = physicalCoordinate;

    case 'xi'
        if ~isfield(solution, 'xi')
            error('extractOWNSCoordinate:MissingXi', ...
                'solution.xi is unavailable.');
        end

        xCoordinate = solution.xi(:);
        physicalCoordinate = xCoordinate;

        if numel(xCoordinate) ~= Nx
            error('extractOWNSCoordinate:XiLengthMismatch', ...
                'solution.xi does not contain Nx entries.');
        end

    otherwise
        error('extractOWNSCoordinate:UnknownMethod', ...
            'Unknown coordinate method "%s".', ...
            coordinateCfg.method);
end

if coordinateCfg.zeroOrigin
    xCoordinate = xCoordinate - xCoordinate(1);
end

if any(~isfinite(xCoordinate))
    error('extractOWNSCoordinate:NonfiniteCoordinate', ...
        'The extracted coordinate contains nonfinite values.');
end

increments = diff(xCoordinate);

if any(increments <= 0)
    badIndex = find(increments <= 0, 1, 'first');

    error('extractOWNSCoordinate:NonIncreasingCoordinate', ...
        ['The extracted first-passage coordinate is not strictly ', ...
         'increasing between stations %d and %d.'], ...
        badIndex, badIndex + 1);
end

coordinateMeta = struct();
coordinateMeta.method = method;
coordinateMeta.referenceIndex = referenceIndex;
coordinateMeta.zeroOrigin = coordinateCfg.zeroOrigin;
coordinateMeta.includeZ = coordinateCfg.includeZ;
coordinateMeta.physicalReferenceCoordinate = physicalCoordinate;
coordinateMeta.minimumIncrement = min(increments);
coordinateMeta.maximumIncrement = max(increments);

end

function [arcLength, physicalCoordinate] = ...
    computeArcLength(solution, referenceIndex, includeZ)

xRef = solution.x(referenceIndex, :).';
yRef = solution.y(referenceIndex, :).';

if includeZ && isfield(solution, 'z') && ...
        ~isempty(solution.z)

    zRef = solution.z(referenceIndex, :).';
else
    zRef = zeros(size(xRef));
end

physicalCoordinate = [xRef, yRef, zRef];

increments = sqrt( ...
    diff(xRef).^2 ...
    + diff(yRef).^2 ...
    + diff(zRef).^2);

arcLength = [0; cumsum(increments)];

end
