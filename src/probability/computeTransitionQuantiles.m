function quantiles = computeTransitionQuantiles(x, F, levels)
%COMPUTETRANSITIONQUANTILES Compute generalized inverse-CDF quantiles.
%
% Quantiles not reached within the computational domain are returned as
% NaN and marked as censored.

x = x(:);
F = F(:);
levels = levels(:);

if numel(x) ~= numel(F)
    error('computeTransitionQuantiles:LengthMismatch', ...
        'x and F must have the same length.');
end

if any(diff(x) <= 0)
    error('computeTransitionQuantiles:InvalidGrid', ...
        'x must be strictly increasing.');
end

if any(diff(F) < -1e-12)
    error('computeTransitionQuantiles:NonmonotoneCdf', ...
        'F must be nondecreasing.');
end

nLevels = numel(levels);

values = nan(nLevels, 1);
isCensored = false(nLevels, 1);

for j = 1:nLevels
    p = levels(j);

    if p < 0 || p > 1
        error('computeTransitionQuantiles:InvalidLevel', ...
            'Quantile levels must lie in [0,1].');
    end

    if p > F(end)
        isCensored(j) = true;
        continue;
    end

    index = find(F >= p, 1, 'first');

    if isempty(index)
        isCensored(j) = true;

    elseif index == 1
        values(j) = x(1);

    elseif F(index) == F(index - 1)
        % A flat CDF segment carries no probability. Return the right
        % endpoint selected by the generalized inverse.
        values(j) = x(index);

    else
        fraction = ...
            (p - F(index - 1)) / ...
            (F(index) - F(index - 1));

        values(j) = ...
            x(index - 1) ...
            + fraction * (x(index) - x(index - 1));
    end
end

quantiles = struct();
quantiles.levels = levels;
quantiles.values = values;
quantiles.isCensored = isCensored;
quantiles.terminalCdf = F(end);

end
