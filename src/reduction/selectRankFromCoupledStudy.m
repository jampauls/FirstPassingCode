function selection = selectRankFromCoupledStudy(study, options)
%SELECTRANKFROMCOUPLEDSTUDY Select the smallest acceptable reduced rank.
%
% A rank is accepted when:
%
%   maximumReferenceDifference
%       + confidenceMultiplier * maximumReferenceDifferenceSE
%       <= cdfTolerance
%
% and, optionally:
%
%   relativeFullResidual <= matrixResidualTolerance.
%
% The reference is the largest rank included in the coupled study.
%
% Inputs:
%   study    output of computeCoupledRankStudy
%   options  structure with optional fields:
%
%       cdfTolerance
%           Absolute tolerance for the complete transition CDF.
%
%       confidenceMultiplier
%           Multiplier applied to the paired standard error.
%
%       matrixResidualTolerance
%           Maximum allowed relative matrix residual. Set to Inf to
%           disable the matrix criterion.
%
%       excludeReferenceRank
%           If true, do not select the largest tested rank unless no
%           reduced rank satisfies the criteria.
%
% Output:
%   selection.selectedRank
%   selection.selectedIndex
%   selection.accepted
%   selection.probabilityUpperBound
%   selection.criteria
%
% MATLAB version: R2020b

if nargin < 2
    options = struct();
end

cdfTolerance = getOption( ...
    options, 'cdfTolerance', 5e-3);

confidenceMultiplier = getOption( ...
    options, 'confidenceMultiplier', 2);

matrixResidualTolerance = getOption( ...
    options, 'matrixResidualTolerance', 1e-2);

excludeReferenceRank = getOption( ...
    options, 'excludeReferenceRank', true);

rankValues = study.rankValues(:);

probabilityUpperBound = ...
    study.maximumReferenceDifference(:) ...
    + confidenceMultiplier ...
    * study.maximumReferenceDifferenceSE(:);

matrixResidual = ...
    study.relativeFullResidual(:);

probabilityCriterion = ...
    probabilityUpperBound <= cdfTolerance;

matrixCriterion = ...
    matrixResidual <= matrixResidualTolerance;

accepted = ...
    probabilityCriterion & matrixCriterion;

candidateIndices = (1:numel(rankValues)).';

if excludeReferenceRank && numel(candidateIndices) > 1
    candidateIndices = candidateIndices(1:end - 1);
end

selectedIndex = [];

for j = candidateIndices.'
    if accepted(j)
        selectedIndex = j;
        break;
    end
end

usedReferenceRank = false;

if isempty(selectedIndex)
    % Fall back to the largest tested rank.
    selectedIndex = numel(rankValues);
    usedReferenceRank = true;
end

selection = struct();

selection.selectedRank = ...
    rankValues(selectedIndex);

selection.selectedIndex = selectedIndex;

selection.accepted = ...
    accepted(selectedIndex);

selection.usedReferenceRank = ...
    usedReferenceRank;

selection.probabilityUpperBound = ...
    probabilityUpperBound;

selection.probabilityCriterion = ...
    probabilityCriterion;

selection.matrixCriterion = ...
    matrixCriterion;

selection.criteria = struct();
selection.criteria.cdfTolerance = cdfTolerance;
selection.criteria.confidenceMultiplier = confidenceMultiplier;
selection.criteria.matrixResidualTolerance = ...
    matrixResidualTolerance;
selection.criteria.excludeReferenceRank = excludeReferenceRank;

fprintf('\nReduced-rank selection:\n');
fprintf('  selected rank                 = %d\n', ...
    selection.selectedRank);
fprintf('  paired CDF difference         = %.6e\n', ...
    study.maximumReferenceDifference(selectedIndex));
fprintf('  paired standard error         = %.6e\n', ...
    study.maximumReferenceDifferenceSE(selectedIndex));
fprintf('  confidence-adjusted CDF error = %.6e\n', ...
    probabilityUpperBound(selectedIndex));
fprintf('  relative matrix residual      = %.6e\n', ...
    matrixResidual(selectedIndex));
fprintf('  CDF tolerance                 = %.6e\n', ...
    cdfTolerance);
fprintf('  matrix residual tolerance     = %.6e\n', ...
    matrixResidualTolerance);

if usedReferenceRank
    fprintf(['  No tested reduced rank satisfied all criteria; ', ...
        'the reference rank was selected.\n']);
elseif selection.accepted
    fprintf('  All configured selection criteria were satisfied.\n');
end

fprintf('\n');

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end
