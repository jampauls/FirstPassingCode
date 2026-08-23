function diagnostics = validateProbabilityCurves(prob, cfg)
%VALIDATEPROBABILITYCURVES Validate survival and first-transition output.

if nargin < 2
    cfg = struct();
end

S = prob.S(:);
F = prob.F(:);
pInterval = prob.pInterval(:);

if numel(S) ~= numel(F)
    error('validateProbabilityCurves:LengthMismatch', ...
        'S and F must have the same length.');
end

if any(~isfinite(S)) || any(~isfinite(F))
    error('validateProbabilityCurves:NonfiniteProbability', ...
        'S or F contains nonfinite values.');
end

tol = 1e-10;

if any(S < -tol) || any(S > 1 + tol)
    error('validateProbabilityCurves:SurvivalBounds', ...
        'Survival probabilities lie outside [0,1].');
end

if any(F < -tol) || any(F > 1 + tol)
    error('validateProbabilityCurves:CdfBounds', ...
        'Transition CDF values lie outside [0,1].');
end

complementarityDefect = max(abs(S + F - 1));

if complementarityDefect > tol
    error('validateProbabilityCurves:Complementarity', ...
        'Maximum S + F - 1 defect is %.3e.', ...
        complementarityDefect);
end

if numel(F) > 1
    minCdfIncrement = min(diff(F));
    maxSurvivalIncrement = max(diff(S));
else
    minCdfIncrement = 0;
    maxSurvivalIncrement = 0;
end

if minCdfIncrement < -tol
    error('validateProbabilityCurves:NonmonotoneCdf', ...
        'The transition CDF decreases by %.3e.', ...
        minCdfIncrement);
end

if maxSurvivalIncrement > tol
    error('validateProbabilityCurves:NonmonotoneSurvival', ...
        'The survival probability increases by %.3e.', ...
        maxSurvivalIncrement);
end

if any(pInterval < -tol)
    error('validateProbabilityCurves:NegativeIntervalProbability', ...
        'An interval probability is negative.');
end

totalProbability = sum(pInterval) + prob.pCensored;
conservationDefect = abs(totalProbability - 1);

if conservationDefect > tol
    error('validateProbabilityCurves:ProbabilityConservation', ...
        'Probability conservation defect is %.3e.', ...
        conservationDefect);
end

diagnostics = struct();
diagnostics.complementarityDefect = complementarityDefect;
diagnostics.minCdfIncrement = minCdfIncrement;
diagnostics.maxSurvivalIncrement = maxSurvivalIncrement;
diagnostics.conservationDefect = conservationDefect;
diagnostics.cfg = cfg;

end
