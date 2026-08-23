function diagnostics = validateRunningMax(m)
%VALIDATERUNNINGMAX Validate directional running maxima.

if ~isnumeric(m) || ndims(m) ~= 2 || isempty(m)
    error('validateRunningMax:InvalidInput', ...
        'Running maxima must be a nonempty numeric matrix.');
end

if any(~isfinite(m(:)))
    error('validateRunningMax:NonfiniteValues', ...
        'Running maxima contain nonfinite values.');
end

minimumValue = min(m(:));

if minimumValue < -1e-12
    error('validateRunningMax:NegativeMaximum', ...
        'Running maxima contain a negative value %.3e.', minimumValue);
end

if size(m, 2) > 1
    increments = diff(m, 1, 2);
    minimumIncrement = min(increments(:));
else
    minimumIncrement = 0;
end

if minimumIncrement < -1e-12
    error('validateRunningMax:Nonmonotone', ...
        'Running maxima decrease by as much as %.3e.', ...
        minimumIncrement);
end

diagnostics = struct();
diagnostics.minimumValue = minimumValue;
diagnostics.minimumIncrement = minimumIncrement;

end