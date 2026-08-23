function replicate = evaluateAngularReplicate( ...
    x, G, Gprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%EVALUATEANGULARREPLICATE Evaluate one angular integration replicate.
%
% This function generates one angular rule, computes directional running
% maxima, and evaluates the resulting first-transition CDF.
%
% Inputs:
%   x               streamwise grid
%   G               normalized energy matrices
%   Gprime          derivative matrices, or []
%   r               stochastic dimension
%   angularCfg      angular-rule configuration
%   maxDetectionCfg maximum-detection configuration
%   numericsCfg     numerical configuration
%
% Output:
%   replicate.U
%   replicate.weights
%   replicate.angularMeta
%   replicate.maxData
%   replicate.prob
%
% MATLAB version: R2020b

cfg = struct();
cfg.angular = angularCfg;
cfg.maxDetection = maxDetectionCfg;
cfg.numerics = numericsCfg;
cfg.gaussianConvention = 'real';

[U, weights, angularMeta] = ...
    generateAngularRule(r, angularCfg);

validateAngularRule(U, weights);

switch lower(maxDetectionCfg.method)

    case 'grid'
        maxData = computeRunningMaxGrid(G, U, cfg);

    case 'hermite'
        if isempty(Gprime)
            error('evaluateAngularReplicate:MissingGprime', ...
                'Hermite maximum detection requires Gprime.');
        end

        maxData = computeRunningMaxHermite( ...
            x, G, Gprime, U, cfg);

    case 'adaptive'
        if isempty(Gprime)
            error('evaluateAngularReplicate:MissingGprime', ...
                'Adaptive maximum detection requires Gprime.');
        end

        maxData = computeRunningMaxAdaptive( ...
            x, G, Gprime, U, cfg);

    otherwise
        error('evaluateAngularReplicate:UnknownMaximumMethod', ...
            'Unknown maximum method "%s".', ...
            maxDetectionCfg.method);
end

validateRunningMax(maxData.m);

prob = computeTransitionCDF( ...
    maxData.m, weights, r, cfg);

validateProbabilityCurves(prob, cfg);

replicate = struct();
replicate.U = U;
replicate.weights = weights;
replicate.angularMeta = angularMeta;
replicate.maxData = maxData;
replicate.prob = prob;

end
