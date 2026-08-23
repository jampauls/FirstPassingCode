function envelope = buildRQMCEnergyEnvelope( ...
    x, A, Aprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%BUILDRQMCENERGYENVELOPE Build unnormalized directional energy envelopes.
%
%   envelope = buildRQMCEnergyEnvelope( ...
%       x, A, Aprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%
% The matrices A are not divided by a transition threshold. For each RQMC
% replicate and direction, this function computes:
%
%   a_k(x_n) = u_k' * A(x_n) * u_k
%
% and
%
%   m_k(x_n) = max_{xi <= x_n} a_k(xi).
%
% Once these quantities are available, many constant thresholds and inlet
% amplitudes can be evaluated without repeating directional maximization.
%
% Inputs:
%   x               Nx-by-1 streamwise coordinate
%   A               cell array of r-by-r energy matrices
%   Aprime          cell array of A'(x), or []
%   r               stochastic dimension
%   angularCfg      RQMC configuration
%   maxDetectionCfg maximum-detection configuration
%   numericsCfg     numerical configuration
%
% Output:
%   envelope.mReplicate       cell{J}, each K-by-Nx
%   envelope.aReplicate       cell{J}, each K-by-Nx
%   envelope.weights          K-by-1
%   envelope.directionMeta    cell{J}
%   envelope.x
%
% MATLAB version: R2020b

x = x(:);
Nx = numel(x);

if numel(A) ~= Nx
    error('buildRQMCEnergyEnvelope:LengthMismatch', ...
        'A must contain one matrix per streamwise station.');
end

if ~isempty(Aprime) && numel(Aprime) ~= Nx
    error('buildRQMCEnergyEnvelope:DerivativeLengthMismatch', ...
        'Aprime must be empty or contain one matrix per station.');
end

if ~strcmpi(angularCfg.method, 'rqmc')
    error('buildRQMCEnergyEnvelope:AngularMethod', ...
        'The current implementation requires angularCfg.method = ''rqmc''.');
end

K = angularCfg.numDirections;
J = angularCfg.numReplicates;

validateattributes(K, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'angularCfg.numDirections');

validateattributes(J, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, 'angularCfg.numReplicates');

weights = ones(K, 1) / K;

mReplicate = cell(J, 1);
aReplicate = cell(J, 1);
directionMeta = cell(J, 1);
maximumInfo = cell(J, 1);

cfgLocal = struct();
cfgLocal.numerics = numericsCfg;
cfgLocal.maxDetection = maxDetectionCfg;
cfgLocal.gaussianConvention = 'real';

fprintf('\nBuilding reusable RQMC energy envelope:\n');
fprintf('  stochastic dimension = %d\n', r);
fprintf('  directions/replicate = %d\n', K);
fprintf('  replicates           = %d\n', J);
fprintf('  maximum method       = %s\n\n', ...
    maxDetectionCfg.method);

for j = 1:J
    fprintf('  envelope replicate %d of %d\n', j, J);

    replicateCfg = angularCfg;
    replicateCfg.replicateIndex = j;

    [U, replicateWeights, angularMeta] = ...
        generateAngularRule(r, replicateCfg);

    validateAngularRule(U, replicateWeights);

    switch lower(maxDetectionCfg.method)

        case 'grid'
            maxData = computeRunningMaxGrid( ...
                A, U, cfgLocal);

        case 'hermite'
            if isempty(Aprime)
                error('buildRQMCEnergyEnvelope:MissingAprime', ...
                    'Hermite maximum detection requires Aprime.');
            end

            maxData = computeRunningMaxHermite( ...
                x, A, Aprime, U, cfgLocal);

        case 'adaptive'
            if isempty(Aprime)
                error('buildRQMCEnergyEnvelope:MissingAprime', ...
                    'Adaptive maximum detection requires Aprime.');
            end

            maxData = computeRunningMaxAdaptive( ...
                x, A, Aprime, U, cfgLocal);

        otherwise
            error('buildRQMCEnergyEnvelope:UnknownMaximumMethod', ...
                'Unknown maximum method "%s".', ...
                maxDetectionCfg.method);
    end

    validateRunningMax(maxData.m);

    mReplicate{j} = maxData.m;
    aReplicate{j} = maxData.a;
    directionMeta{j} = angularMeta;
    maximumInfo{j} = maxData.info;
end

envelope = struct();

envelope.x = x;
envelope.r = r;
envelope.numStations = Nx;
envelope.numDirections = K;
envelope.numReplicates = J;

envelope.weights = weights;
envelope.mReplicate = mReplicate;
envelope.aReplicate = aReplicate;

envelope.directionMeta = directionMeta;
envelope.maximumInfo = maximumInfo;

envelope.angularCfg = angularCfg;
envelope.maxDetectionCfg = maxDetectionCfg;
envelope.numericsCfg = numericsCfg;

end
