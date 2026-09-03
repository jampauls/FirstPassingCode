function indices = selectLocalStencil( ...
    numStations, centerIndex, stencilSize)
%SELECTLOCALSTENCIL Select a centered or one-sided local stencil.
%
% MATLAB version: R2020b

validateattributes(numStations, {'numeric'}, ...
    {'scalar', 'integer', 'positive'});

validateattributes(centerIndex, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, '<=', numStations});

validateattributes(stencilSize, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2, '<=', numStations});

leftCount = floor((stencilSize - 1) / 2);

firstIndex = centerIndex - leftCount;
lastIndex = firstIndex + stencilSize - 1;

if firstIndex < 1
    firstIndex = 1;
    lastIndex = stencilSize;
elseif lastIndex > numStations
    lastIndex = numStations;
    firstIndex = numStations - stencilSize + 1;
end

indices = firstIndex:lastIndex;

end