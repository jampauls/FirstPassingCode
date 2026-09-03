function [tailProbability, info] = generalizedQuadraticFormTail( ...
    eigenvalues, level, gaussianConvention, options)
%GENERALIZEDQUADRATICFORMTAIL Upper tail of a Gaussian quadratic form.
%
% Real Gaussian case:
%
%   Q = sum_j lambda_j * Z_j^2,
%   Z_j ~ N(0,1).
%
% Proper-complex Gaussian case:
%
%   Q = sum_j lambda_j * |Z_j|^2,
%   Z_j ~ CN(0,1),
%   |Z_j|^2 ~ Exponential(1).
%
% The function calculates:
%
%   tailProbability = P(Q >= level)
%
% using Gil--Pelaez characteristic-function inversion. Equal-eigenvalue
% and one-active-eigenvalue cases are evaluated analytically.
%
% Positive-semidefinite, negative-semidefinite, and indefinite quadratic
% forms are supported. Characteristic-function inversion is used for
% general signed eigenvalue families.

%
% Options:
%   absoluteTolerance          default 1e-10
%   relativeTolerance          default 1e-8
%   eigenvalueTolerance        default 1e-13
%   equalEigenvalueTolerance   default 1e-12
%   maxIntervalCount           default 10000
%   warnOnLargeError           default true
%
% MATLAB version: R2020b

if nargin < 3 || isempty(gaussianConvention)
    gaussianConvention = 'real';
end

if nargin < 4 || isempty(options)
    options = struct();
end

absoluteTolerance = getOption( ...
    options, 'absoluteTolerance', 1e-10);

relativeTolerance = getOption( ...
    options, 'relativeTolerance', 1e-8);

eigenvalueTolerance = getOption( ...
    options, 'eigenvalueTolerance', 1e-13);

equalEigenvalueTolerance = getOption( ...
    options, 'equalEigenvalueTolerance', 1e-12);

maxIntervalCount = getOption( ...
    options, 'maxIntervalCount', 10000);

warnOnLargeError = getOption( ...
    options, 'warnOnLargeError', true);

validateattributes(level, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, ...
    mfilename, 'level');

lambda = real(eigenvalues(:));

if any(~isfinite(lambda))
    error('generalizedQuadraticFormTail:NonfiniteEigenvalue', ...
        'All eigenvalues must be finite.');
end

% -------------------------------------------------------------------------
% Remove roundoff-scale eigenvalues of either sign
% -------------------------------------------------------------------------

lambdaScale = max(max(abs(lambda)), 1);

smallEigenvalue = ...
    abs(lambda) <= eigenvalueTolerance * lambdaScale;

lambda(smallEigenvalue) = 0;

active = lambda ~= 0;
lambda = lambda(active);

numActive = numel(lambda);

positiveLambda = lambda(lambda > 0);
negativeLambda = lambda(lambda < 0);

info = struct();
info.gaussianConvention = char(gaussianConvention);
info.numActiveEigenvalues = numActive;
info.numPositiveEigenvalues = numel(positiveLambda);
info.numNegativeEigenvalues = numel(negativeLambda);
info.isIndefinite = ...
    ~isempty(positiveLambda) && ~isempty(negativeLambda);
info.usedAnalyticalFormula = false;
info.integrationErrorEstimate = 0;
info.scalingValue = NaN;

% -------------------------------------------------------------------------
% Degenerate and support-boundary cases
% -------------------------------------------------------------------------

if numActive == 0
    tailProbability = double(level <= 0);
    info.method = 'degenerate';
    return;
end

% If all active eigenvalues are nonnegative, Q >= 0.
if isempty(negativeLambda) && level <= 0
    tailProbability = 1;
    info.method = 'nonnegativeSupport';
    return;
end

% If all active eigenvalues are nonpositive, Q <= 0.
if isempty(positiveLambda) && level > 0
    tailProbability = 0;
    info.method = 'nonpositiveSupport';
    return;
end
% -------------------------------------------------------------------------
% Scale the eigenvalues and level
% -------------------------------------------------------------------------

scale = max(abs(lambda));

lambdaScaled = lambda / scale;
levelScaled = level / scale;

info.scalingValue = scale;

% -------------------------------------------------------------------------
% Analytical one-eigenvalue case
% -------------------------------------------------------------------------

if numActive == 1
    lambdaOne = lambda(1);

    switch lower(char(gaussianConvention))

        case 'real'
            if lambdaOne > 0
                if level <= 0
                    tailProbability = 1;
                else
                    tailProbability = gammainc( ...
                        level / (2 * lambdaOne), ...
                        1 / 2, ...
                        'upper');
                end
            else
                if level > 0
                    tailProbability = 0;
                else
                    tailProbability = gammainc( ...
                        (-level) / (2 * abs(lambdaOne)), ...
                        1 / 2, ...
                        'lower');
                end
            end

            info.method = 'singleEigenvalueChiSquare';

        case 'propercomplex'
            if lambdaOne > 0
                if level <= 0
                    tailProbability = 1;
                else
                    tailProbability = ...
                        exp(-level / lambdaOne);
                end
            else
                if level > 0
                    tailProbability = 0;
                else
                    tailProbability = ...
                        1 - exp(level / abs(lambdaOne));
                end
            end

            info.method = 'singleEigenvalueExponential';

        otherwise
            error('generalizedQuadraticFormTail:UnknownConvention', ...
                'Unknown Gaussian convention "%s".', ...
                gaussianConvention);
    end

    info.usedAnalyticalFormula = true;
    tailProbability = min(max(real(tailProbability), 0), 1);
    return;
end

% -------------------------------------------------------------------------
% Equal-eigenvalue analytical case
% -------------------------------------------------------------------------

commonEigenvalue = mean(lambda);

relativeSpread = ...
    max(abs(lambda - commonEigenvalue)) ...
    / max(abs(commonEigenvalue), eps);

if relativeSpread <= equalEigenvalueTolerance
    switch lower(char(gaussianConvention))

        case 'real'
            if commonEigenvalue > 0
                if level <= 0
                    tailProbability = 1;
                else
                    tailProbability = gammainc( ...
                        level / (2 * commonEigenvalue), ...
                        numActive / 2, ...
                        'upper');
                end
            else
                if level > 0
                    tailProbability = 0;
                else
                    tailProbability = gammainc( ...
                        (-level) ...
                            / (2 * abs(commonEigenvalue)), ...
                        numActive / 2, ...
                        'lower');
                end
            end

            info.method = 'equalEigenvalueChiSquare';

        case 'propercomplex'
            if commonEigenvalue > 0
                if level <= 0
                    tailProbability = 1;
                else
                    tailProbability = gammainc( ...
                        level / commonEigenvalue, ...
                        numActive, ...
                        'upper');
                end
            else
                if level > 0
                    tailProbability = 0;
                else
                    tailProbability = gammainc( ...
                        (-level) / abs(commonEigenvalue), ...
                        numActive, ...
                        'lower');
                end
            end

            info.method = 'equalEigenvalueGamma';

        otherwise
            error('generalizedQuadraticFormTail:UnknownConvention', ...
                'Unknown Gaussian convention "%s".', ...
                gaussianConvention);
    end

    info.usedAnalyticalFormula = true;
    tailProbability = min(max(real(tailProbability), 0), 1);
    return;
end


% -------------------------------------------------------------------------
% Characteristic-function inversion
%
% Gil--Pelaez:
%
%   P(Q >= q)
%       = 1/2
%         + (1/pi) integral_0^infinity
%             Im[exp(-itq) phi_Q(t)] / t dt.
% -------------------------------------------------------------------------

switch lower(char(gaussianConvention))

    case 'real'
        integrand = @(t) realGaussianIntegrand( ...
            t, lambdaScaled, levelScaled);

        info.method = 'characteristicFunctionReal';

    case 'propercomplex'
        integrand = @(t) properComplexIntegrand( ...
            t, lambdaScaled, levelScaled);

        info.method = 'characteristicFunctionProperComplex';

    otherwise
        error('generalizedQuadraticFormTail:UnknownConvention', ...
            'Unknown Gaussian convention "%s".', ...
            gaussianConvention);
end

[integralValue, integrationError] = quadgk( ...
    integrand, ...
    0, ...
    Inf, ...
    'AbsTol', absoluteTolerance, ...
    'RelTol', relativeTolerance, ...
    'MaxIntervalCount', maxIntervalCount);

tailProbability = ...
    0.5 + integralValue / pi;

tailProbability = ...
    min(max(real(tailProbability), 0), 1);

info.integrationErrorEstimate = ...
    integrationError / pi;

if warnOnLargeError
    permittedError = max( ...
        10 * absoluteTolerance, ...
        10 * relativeTolerance ...
            * max(tailProbability, eps));

    if info.integrationErrorEstimate > permittedError
        warning( ...
            'generalizedQuadraticFormTail:LargeIntegrationError', ...
            ['Estimated probability integration error %.3e exceeds ', ...
             'the configured diagnostic tolerance %.3e.'], ...
            info.integrationErrorEstimate, permittedError);
    end
end

end

% =========================================================================
function value = realGaussianIntegrand(t, lambda, level)
%REALGAUSSIANINTEGRAND Inversion integrand for weighted chi-square terms.

originalSize = size(t);
tColumn = t(:);

arguments = ...
    2 * tColumn * lambda.';

% Characteristic function:
%
%   phi(t) = prod_j (1 - 2i*t*lambda_j)^(-1/2).
%
% Evaluate through logarithmic amplitude and phase.
logAmplitude = ...
    -0.25 * sum(log1p(arguments.^2), 2);

phase = ...
    0.5 * sum(atan(arguments), 2) ...
    - tColumn * level;

valueColumn = zeros(size(tColumn));

nonzero = tColumn ~= 0;

valueColumn(nonzero) = ...
    exp(logAmplitude(nonzero)) ...
    .* sin(phase(nonzero)) ...
    ./ tColumn(nonzero);

% Limit as t -> 0.
valueColumn(~nonzero) = ...
    sum(lambda) - level;

value = reshape(valueColumn, originalSize);

end

% =========================================================================
function value = properComplexIntegrand(t, lambda, level)
%PROPERCOMPLEXINTEGRAND Inversion integrand for weighted exponentials.

originalSize = size(t);
tColumn = t(:);

arguments = ...
    tColumn * lambda.';

% Characteristic function:
%
%   phi(t) = prod_j (1 - i*t*lambda_j)^(-1).
logAmplitude = ...
    -0.5 * sum(log1p(arguments.^2), 2);

phase = ...
    sum(atan(arguments), 2) ...
    - tColumn * level;

valueColumn = zeros(size(tColumn));

nonzero = tColumn ~= 0;

valueColumn(nonzero) = ...
    exp(logAmplitude(nonzero)) ...
    .* sin(phase(nonzero)) ...
    ./ tColumn(nonzero);

valueColumn(~nonzero) = ...
    sum(lambda) - level;

value = reshape(valueColumn, originalSize);

end

% =========================================================================
function value = getOption(options, fieldName, defaultValue)

if isfield(options, fieldName) && ...
        ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end

end