function result = computeRQMCCDFConverged( ...
    x, G, Gprime, r, angularCfg, maxDetectionCfg, numericsCfg)
%COMPUTERQMCCDFCONVERGED Automatically refine RQMC angular resolution.
%
% The number of points per replicate is doubled until:
%
%   max_x |F_K(x) - F_{K/2}(x)| <= cdfTolerance
%
% and
%
%   max_x SE[F_K(x)] <= standardErrorTolerance.
%
% MATLAB version: R2020b

K = angularCfg.initialNumDirections;
Kmax = angularCfg.maxNumDirections;

cdfTolerance = angularCfg.cdfTolerance;
standardErrorTolerance = ...
    angularCfg.standardErrorTolerance;

validateattributes(K, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'initialNumDirections');

validateattributes(Kmax, {'numeric'}, ...
    {'scalar', 'integer', '>=', K}, ...
    mfilename, 'maxNumDirections');

history = struct();
history.K = [];
history.terminalProbability = [];
history.maxCdfChange = [];
history.maxStandardError = [];

previous = [];
converged = false;
iteration = 0;

while K <= Kmax
    iteration = iteration + 1;

    currentCfg = angularCfg;
    currentCfg.method = 'rqmc';
    currentCfg.numDirections = K;

    fprintf('RQMC angular calculation: K = %d, J = %d\n', ...
        K, currentCfg.numReplicates);

    current = computeRQMCCDF( ...
        x, G, Gprime, r, ...
        currentCfg, maxDetectionCfg, numericsCfg);

    maxStandardError = max(current.standardError);

    if isempty(previous)
        maxCdfChange = NaN;
    else
        maxCdfChange = max(abs(current.F - previous.F));
    end

    history.K(iteration, 1) = K;
    history.terminalProbability(iteration, 1) = current.F(end);
    history.maxCdfChange(iteration, 1) = maxCdfChange;
    history.maxStandardError(iteration, 1) = maxStandardError;

    fprintf('  F(xmax)       = %.12e\n', current.F(end));
    fprintf('  terminal SE   = %.12e\n', ...
        current.terminalStandardError);
    fprintf('  max CDF SE    = %.12e\n', maxStandardError);

    if ~isnan(maxCdfChange)
        fprintf('  max CDF change = %.12e\n', maxCdfChange);
    end

    if ~isempty(previous) && ...
            maxCdfChange <= cdfTolerance && ...
            maxStandardError <= standardErrorTolerance

        converged = true;
        break;
    end

    previous = current;
    K = 2 * K;
end

if ~converged
    warning('computeRQMCCDFConverged:NotConverged', ...
        ['RQMC angular calculation did not satisfy the requested ', ...
         'tolerances before reaching Kmax.']);
end

result = current;
result.converged = converged;
result.convergenceHistory = history;
result.finalNumDirections = currentCfg.numDirections;

end
