function [eThresh, info] = ...
    buildEnsembleTransitionThreshold(x, meanEnergyTotal, varianceEnergyTotal, thresholdCfg)
%BUILDENSEMBLETRANSITIONTHRESHOLD Threshold for total ensemble energy.
%
%   [eThresh, info] = buildEnsembleTransitionThreshold( ...
%       x, meanEnergyTotal, varianceEnergyTotal, thresholdCfg)
%
% meanEnergyTotal(x) and varianceEnergyTotal(x) are the mean and variance
% of the energy summed independently over every (omega,beta) mode.
%
% Supported cfg.threshold.method values:
%
%   'specifiedScalar'    : constant, user-supplied absolute threshold
%   'initialMeanMultiple': multiplier * mean total energy at the inlet
%   'auto'               : autoMaxMeanFactor * max mean total energy
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

meanEnergyTotal = meanEnergyTotal(:);
varianceEnergyTotal = varianceEnergyTotal(:);

if numel(meanEnergyTotal) ~= Nx || numel(varianceEnergyTotal) ~= Nx
    error('buildEnsembleTransitionThreshold:LengthMismatch', ...
        'meanEnergyTotal and varianceEnergyTotal must contain Nx entries.');
end

if nargin < 4 || isempty(thresholdCfg)
    thresholdCfg = struct();
end

if ~isfield(thresholdCfg, 'method') || isempty(thresholdCfg.method)
    thresholdCfg.method = 'auto';
end

if ~isfield(thresholdCfg, 'warnIfProvisional')
    thresholdCfg.warnIfProvisional = true;
end

method = lower(char(thresholdCfg.method));
isProvisional = false;

switch method

    case 'specifiedscalar'
        if ~isfield(thresholdCfg, 'value') || isempty(thresholdCfg.value)
            error('buildEnsembleTransitionThreshold:MissingScalar', ...
                'thresholdCfg.value is required.');
        end

        validateattributes(thresholdCfg.value, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});

        eThresh = repmat(thresholdCfg.value, Nx, 1);

    case 'initialmeanmultiple'
        multiplier = getOption(thresholdCfg, 'multiplier', 10);

        thresholdValue = multiplier * meanEnergyTotal(1);

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    case 'auto'
        factor = getOption(thresholdCfg, 'autoMaxMeanFactor', 1);

        thresholdValue = factor * max(meanEnergyTotal);

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    otherwise
        error('buildEnsembleTransitionThreshold:UnknownMethod', ...
            ['Unknown ensemble threshold method "%s". Supported ', ...
             'methods are specifiedScalar, initialMeanMultiple, ', ...
             'and auto.'], thresholdCfg.method);
end

if any(~isfinite(eThresh)) || any(eThresh <= 0)
    error('buildEnsembleTransitionThreshold:InvalidThreshold', ...
        'The generated threshold must be finite and positive.');
end

if isProvisional && thresholdCfg.warnIfProvisional
    warning('buildEnsembleTransitionThreshold:ProvisionalThreshold', ...
        ['Using provisional data-derived ensemble transition threshold ', ...
         'method "%s". This threshold is suitable for numerical ', ...
         'development but is not a validated physical criterion.'], ...
        thresholdCfg.method);
end

info = struct();
info.method = method;
info.isProvisional = isProvisional;
info.meanEnergyTotal = meanEnergyTotal;
info.varianceEnergyTotal = varianceEnergyTotal;
info.minimumThreshold = min(eThresh);
info.maximumThreshold = max(eThresh);
info.initialMeanEnergyTotal = meanEnergyTotal(1);
info.maximumMeanEnergyTotal = max(meanEnergyTotal);

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end
