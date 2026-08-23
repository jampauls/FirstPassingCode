function [thresholdValues, info] = ...
    buildThresholdSweepValues(A, sweepCfg)
%BUILDTHRESHOLDSWEEPVALUES Construct constant threshold values.
%
% MATLAB version: R2020b

Nx = numel(A);

meanEnergy = zeros(Nx, 1);

for n = 1:Nx
    meanEnergy(n) = trace(A{n});
end

if isfield(sweepCfg, 'values') && ...
        ~isempty(sweepCfg.values)

    thresholdValues = sweepCfg.values(:);

    referenceValue = NaN;
    method = 'specifiedValues';

else
    method = lower(char( ...
        sweepCfg.referenceMethod));

    switch method

        case 'maximummeanenergy'
            referenceValue = max(meanEnergy);

        case 'initialmeanenergy'
            referenceValue = meanEnergy(1);

        otherwise
            error('buildThresholdSweepValues:UnknownMethod', ...
                'Unknown referenceMethod "%s".', ...
                sweepCfg.referenceMethod);
    end

    factors = sweepCfg.factors(:);

    thresholdValues = ...
        referenceValue * factors;
end

if any(~isfinite(thresholdValues)) || ...
        any(thresholdValues <= 0)
    error('buildThresholdSweepValues:InvalidThresholds', ...
        'All threshold values must be finite and positive.');
end

thresholdValues = sort(unique(thresholdValues));

info = struct();
info.method = method;
info.referenceValue = referenceValue;
info.meanEnergy = meanEnergy;
info.minimumThreshold = min(thresholdValues);
info.maximumThreshold = max(thresholdValues);

end
