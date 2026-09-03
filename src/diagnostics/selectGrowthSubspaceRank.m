function selection = selectGrowthSubspaceRank(report, options)
%SELECTGROWTHSUBSPACERANK Select a candidate fixed growth-basis rank.
%
% Criteria:
%   aggregateCapturedFraction >= minimumAggregateCapture
%   weightedPositiveTraceCapture >= minimumTraceCapture
%   maximumRelativeComplementGrowth <= maximumComplementGrowth
%   maximumRelativeCrossCoupling <= maximumCrossCoupling
%
% MATLAB version: R2020b

if nargin < 2
    options = struct();
end

minimumAggregateCapture = getOption( ...
    options, 'minimumAggregateCapture', 0.999);

minimumTraceCapture = getOption( ...
    options, 'minimumTraceCapture', 0.999);

maximumComplementGrowth = getOption( ...
    options, 'maximumComplementGrowth', 1e-3);

maximumCrossCoupling = getOption( ...
    options, 'maximumCrossCoupling', 1e-2);

accepted = ...
    report.aggregateCapturedFraction >= minimumAggregateCapture ...
    & report.weightedPositiveTraceCapture >= minimumTraceCapture ...
    & report.maximumRelativeComplementGrowth <= ...
        maximumComplementGrowth ...
    & report.maximumRelativeCrossCoupling <= ...
        maximumCrossCoupling;

index = find(accepted, 1, 'first');

selection = struct();
selection.acceptedByRank = accepted;

if isempty(index)
    selection.found = false;
    selection.selectedRank = [];
    selection.selectedIndex = [];
else
    selection.found = true;
    selection.selectedRank = ...
        report.candidateRanks(index);
    selection.selectedIndex = index;
end

selection.criteria = struct();
selection.criteria.minimumAggregateCapture = ...
    minimumAggregateCapture;
selection.criteria.minimumTraceCapture = ...
    minimumTraceCapture;
selection.criteria.maximumComplementGrowth = ...
    maximumComplementGrowth;
selection.criteria.maximumCrossCoupling = ...
    maximumCrossCoupling;

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end
