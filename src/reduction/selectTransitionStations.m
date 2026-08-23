function [stationIndices, info] = ...
    selectTransitionStations(G, x, reductionCfg)
%SELECTTRANSITIONSTATIONS Select stations for eigenspace snapshots.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(G) ~= Nx
    error('selectTransitionStations:LengthMismatch', ...
        'G must contain one matrix per station.');
end

selection = 'recordsanduniform';

if isfield(reductionCfg, 'stationSelection')
    selection = lower(char( ...
        reductionCfg.stationSelection));
end

numUniform = 10;

if isfield(reductionCfg, 'numUniformStations')
    numUniform = reductionCfg.numUniformStations;
end

recordRelativeTolerance = 1e-6;

if isfield(reductionCfg, 'recordRelativeTolerance')
    recordRelativeTolerance = ...
        reductionCfg.recordRelativeTolerance;
end

lambdaMax = zeros(Nx, 1);

for n = 1:Nx
    lambdaMax(n) = max(eig( ...
        0.5 * (G{n} + G{n}.')));
end

recordIndices = [];
currentRecord = -inf;

for n = 1:Nx
    comparisonScale = max(abs(currentRecord), 1);

    if isempty(recordIndices) || ...
            lambdaMax(n) > currentRecord ...
            + recordRelativeTolerance * comparisonScale

        recordIndices(end + 1, 1) = n; %#ok<AGROW>
        currentRecord = lambdaMax(n);
    end
end

numUniform = min(max(round(numUniform), 2), Nx);

uniformIndices = unique(round( ...
    linspace(1, Nx, numUniform))).';

switch selection

    case 'records'
        stationIndices = recordIndices;

    case 'uniform'
        stationIndices = uniformIndices;

    case 'all'
        stationIndices = (1:Nx).';

    case 'recordsanduniform'
        stationIndices = unique([ ...
            recordIndices; ...
            uniformIndices; ...
            1; ...
            Nx]);

    otherwise
        error('selectTransitionStations:UnknownSelection', ...
            'Unknown stationSelection "%s".', ...
            reductionCfg.stationSelection);
end

stationIndices = sort(stationIndices);

info = struct();
info.selection = selection;
info.stationIndices = stationIndices;
info.stationCoordinates = x(stationIndices);
info.recordIndices = recordIndices;
info.recordCoordinates = x(recordIndices);
info.uniformIndices = uniformIndices;
info.lambdaMax = lambdaMax;

end