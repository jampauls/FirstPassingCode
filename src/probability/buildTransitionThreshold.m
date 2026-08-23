function [eThresh, info] = ...
    buildTransitionThreshold(x, A, thresholdCfg)
%BUILDTRANSITIONTHRESHOLD Construct a configured transition threshold.
%
% This function provides numerical defaults for development and testing.
% Data-derived thresholds are provisional and must eventually be replaced
% or calibrated using DNS or another physical transition criterion.
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(A) ~= Nx
    error('buildTransitionThreshold:LengthMismatch', ...
        'A must contain one matrix per streamwise station.');
end

if nargin < 3 || isempty(thresholdCfg)
    thresholdCfg = struct();
end

if ~isfield(thresholdCfg, 'method') || ...
        isempty(thresholdCfg.method)
    thresholdCfg.method = 'auto';
end

if ~isfield(thresholdCfg, 'warnIfProvisional')
    thresholdCfg.warnIfProvisional = true;
end

meanEnergy = zeros(Nx, 1);
energyVariance = zeros(Nx, 1);

for n = 1:Nx
    An = A{n};

    meanEnergy(n) = trace(An);

    % Real Gaussian quadratic form:
    %
    %   Var(w'Aw) = 2*trace(A^2).
    energyVariance(n) = ...
        2 * trace(An * An);
end

method = lower(char(thresholdCfg.method));
isProvisional = false;

switch method

    case 'specifiedscalar'
        if ~isfield(thresholdCfg, 'value') || ...
                isempty(thresholdCfg.value)
            error('buildTransitionThreshold:MissingScalar', ...
                'thresholdCfg.value is required.');
        end

        validateattributes(thresholdCfg.value, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});

        eThresh = repmat(thresholdCfg.value, Nx, 1);

    case 'specifiedvector'
        if ~isfield(thresholdCfg, 'vector') || ...
                isempty(thresholdCfg.vector)
            error('buildTransitionThreshold:MissingVector', ...
                'thresholdCfg.vector is required.');
        end

        eThresh = thresholdCfg.vector(:);

        if numel(eThresh) ~= Nx
            error('buildTransitionThreshold:VectorLengthMismatch', ...
                'The threshold vector must contain Nx entries.');
        end

    case 'initialmeanmultiple'
        multiplier = getOption( ...
            thresholdCfg, 'multiplier', 10);

        thresholdValue = ...
            multiplier * meanEnergy(1);

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    case 'localmeanmultiple'
        multiplier = getOption( ...
            thresholdCfg, 'multiplier', 10);

        eThresh = multiplier * meanEnergy;
        isProvisional = true;

    case 'maxmeanfraction'
        fraction = getOption( ...
            thresholdCfg, 'fraction', 1);

        thresholdValue = ...
            fraction * max(meanEnergy);

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    case 'initialquantileapprox'
        exceedanceProbability = getOption( ...
            thresholdCfg, ...
            'exceedanceProbability', ...
            1e-3);

        validateattributes(exceedanceProbability, {'numeric'}, ...
            {'scalar', 'real', 'positive', '<', 1});

        mu = meanEnergy(1);
        variance = energyVariance(1);

        if mu <= 0 || variance <= 0
            error('buildTransitionThreshold:DegenerateInitialEnergy', ...
                ['The inlet energy distribution is degenerate and ', ...
                 'cannot be moment-matched to a gamma distribution.']);
        end

        % Satterthwaite/gamma moment matching:
        %
        %   shape = mu^2 / variance
        %   scale = variance / mu
        shape = mu^2 / variance;
        scale = variance / mu;

        lowerProbability = 1 - exceedanceProbability;

        thresholdValue = ...
            scale * gammaincinv(lowerProbability, shape, 'lower');

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    case 'auto'
        % Development default:
        %
        % Set one constant threshold from the maximum mean energy over
        % the complete OWNS domain. This usually gives a numerically
        % nontrivial transition probability while remaining independent
        % of angular sampling.
        factor = getOption( ...
            thresholdCfg, ...
            'autoMaxMeanFactor', ...
            1);

        thresholdValue = ...
            factor * max(meanEnergy);

        eThresh = repmat(thresholdValue, Nx, 1);
        isProvisional = true;

    otherwise
        error('buildTransitionThreshold:UnknownMethod', ...
            'Unknown threshold method "%s".', ...
            thresholdCfg.method);
end

if any(~isfinite(eThresh)) || any(eThresh <= 0)
    error('buildTransitionThreshold:InvalidThreshold', ...
        'The generated threshold must be finite and positive.');
end

if isProvisional && thresholdCfg.warnIfProvisional
    warning('buildTransitionThreshold:ProvisionalThreshold', ...
        ['Using provisional data-derived transition threshold method ', ...
         '"%s". This threshold is suitable for numerical development ', ...
         'but is not a validated physical transition criterion.'], ...
        thresholdCfg.method);
end

info = struct();
info.method = method;
info.isProvisional = isProvisional;
info.meanEnergy = meanEnergy;
info.energyVariance = energyVariance;
info.minimumThreshold = min(eThresh);
info.maximumThreshold = max(eThresh);
info.initialMeanEnergy = meanEnergy(1);
info.maximumMeanEnergy = max(meanEnergy);
info.maximumMeanEnergyIndex = ...
    find(meanEnergy == max(meanEnergy), 1, 'first');
info.maximumMeanEnergyCoordinate = ...
    x(info.maximumMeanEnergyIndex);

end

function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end