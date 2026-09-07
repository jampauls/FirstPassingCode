function solutionOut = downsampleOWNSStations(solution, cfg)
%DOWNSAMPLEOWNSSTATIONS Reduce the number of retained streamwise stations.
%
%   solutionOut = downsampleOWNSStations(solution, cfg)
%
% The OWNS march itself needs on the order of 10,000 streamwise stations
% to resolve its dynamics, but far fewer are needed to sample the
% propagated stochastic structure for first-transition prediction. This
% selects a coarser, still representative subset of stations from the
% raw solution before buildAFromOWNS is ever called, so the expensive
% per-station energy-matrix construction -- and the resulting cache file
% -- only cover the retained stations.
%
% Optional cfg fields:
%   cfg.owns.numCacheStations
%       Target number of retained stations. Empty/absent/>=Nx keeps every
%       station (no downsampling). Tunable per problem.
%
%   cfg.owns.cacheXCutoffCoefficient
%       If set, only stations with
%           (x - x(1)) <= cacheXCutoffCoefficient / sqrt(|omega|)
%       are eligible before the numCacheStations subset is chosen. This
%       concentrates resolution near the inlet for high-frequency modes,
%       whose relevant dynamics occur over that same shrinking length
%       scale, instead of spreading the same station budget across the
%       full, mostly-irrelevant downstream domain. Empty/absent uses the
%       full streamwise domain (no cutoff).
%
% The first and last eligible stations are always retained so the
% downsampled grid still spans its full eligible range.
%
% MATLAB version: R2020b

numStations = [];
if isfield(cfg, 'owns') && isfield(cfg.owns, 'numCacheStations')
    numStations = cfg.owns.numCacheStations;
end

if ~isempty(numStations) && numStations < 2
    error('downsampleOWNSStations:TooFewStations', ...
        'cfg.owns.numCacheStations must be at least 2.');
end

Nx = size(solution.q, 2);
eligible = (1:Nx).';

cutoffCoefficient = [];
if isfield(cfg, 'owns') && isfield(cfg.owns, 'cacheXCutoffCoefficient')
    cutoffCoefficient = cfg.owns.cacheXCutoffCoefficient;
end

% The cutoff is applied first, regardless of numCacheStations, since it
% restricts which stations are even eligible before any station-count
% budget is spent on them.
if ~isempty(cutoffCoefficient)
    if ~isfield(solution, 'w')
        error('downsampleOWNSStations:MissingFrequency', ...
            'solution.w is required to apply a 1/sqrt(omega) cutoff.');
    end

    xWall = solution.x(1, :).';
    xRelative = xWall - xWall(1);

    xCutoff = cutoffCoefficient / sqrt(abs(solution.w));

    withinCutoff = eligible(xRelative <= xCutoff);

    if numel(withinCutoff) >= 2
        eligible = withinCutoff;
    end
end

if isempty(numStations)
    keepCount = numel(eligible);
else
    keepCount = min(numStations, numel(eligible));
end

selected = eligible(round(linspace(1, numel(eligible), keepCount)));
selected = unique(selected, 'stable');

if isequal(selected, (1:Nx).')
    % Every raw station is retained: no cutoff and no downsampling.
    solutionOut = solution;
    return;
end

solutionOut = solution;
solutionOut.q = solution.q(:, selected, :);
solutionOut.gram_W = solution.gram_W(:, selected);
solutionOut.x = solution.x(:, selected);
solutionOut.y = solution.y(:, selected);

if isfield(solution, 'z') && ~isempty(solution.z)
    solutionOut.z = solution.z(:, selected);
end

if isfield(solution, 'xi') && numel(solution.xi) == Nx
    solutionOut.xi = solution.xi(selected);
end

if isfield(solution, 'fp') && ~isempty(solution.fp) && ...
        size(solution.fp, 2) == Nx
    solutionOut.fp = solution.fp(:, selected);
end

end

