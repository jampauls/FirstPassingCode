function [maximumValue, maximumX, info] = ...
    cubicHermiteIntervalMaximum( ...
    xLeft, xRight, valueLeft, valueRight, ...
    derivativeLeft, derivativeRight, options)
%CUBICHERMITEINTERVALMAXIMUM Maximum of a cubic Hermite interpolant.
%
%   [maximumValue, maximumX, info] = ...
%       cubicHermiteIntervalMaximum( ...
%           xLeft, xRight, valueLeft, valueRight, ...
%           derivativeLeft, derivativeRight, options)
%
% All value and derivative inputs may be vectors of equal length. The
% function computes one cubic Hermite maximum per vector entry.
%
% Outputs:
%   maximumValue  maximum interpolated value on [xLeft,xRight]
%   maximumX      estimated maximizing streamwise coordinate
%   info          stationary-point and clipping diagnostics
%
% MATLAB version: R2020b

if nargin < 7
    options = struct();
end

if ~isfield(options, 'polynomialTolerance')
    options.polynomialTolerance = 1e-13;
end

if ~isfield(options, 'rootTolerance')
    options.rootTolerance = 1e-12;
end

if ~isfield(options, 'clipNegativeMaximum')
    options.clipNegativeMaximum = true;
end

validateattributes(xLeft, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, ...
    mfilename, 'xLeft');

validateattributes(xRight, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', xLeft}, ...
    mfilename, 'xRight');

valueLeft = valueLeft(:);
valueRight = valueRight(:);
derivativeLeft = derivativeLeft(:);
derivativeRight = derivativeRight(:);

K = numel(valueLeft);

if numel(valueRight) ~= K || ...
        numel(derivativeLeft) ~= K || ...
        numel(derivativeRight) ~= K
    error('cubicHermiteIntervalMaximum:LengthMismatch', ...
        'All value and derivative inputs must have equal lengths.');
end

if any(~isfinite(valueLeft)) || ...
        any(~isfinite(valueRight)) || ...
        any(~isfinite(derivativeLeft)) || ...
        any(~isfinite(derivativeRight))
    error('cubicHermiteIntervalMaximum:NonfiniteInput', ...
        'Hermite data must be finite.');
end

h = xRight - xLeft;

c0 = valueLeft;
c1 = h * derivativeLeft;
c2 = 3 * (valueRight - valueLeft) ...
    - h * (2 * derivativeLeft + derivativeRight);
c3 = 2 * (valueLeft - valueRight) ...
    + h * (derivativeLeft + derivativeRight);

maximumValue = zeros(K, 1);
maximumX = zeros(K, 1);

numInteriorCandidates = zeros(K, 1);
numNegativeClipped = 0;

for k = 1:K
    coefficients = [c0(k), c1(k), c2(k), c3(k)];

    coefficientScale = max([abs(coefficients), 1]);
    polynomialTolerance = ...
        options.polynomialTolerance * coefficientScale;

    candidates = zeros(4, 1);
    candidates(1:2) = [0; 1];
    numCandidates = 2;

    % Derivative:
    %
    %   3*c3*s^2 + 2*c2*s + c1 = 0.
    A = 3 * c3(k);
    B = 2 * c2(k);
    C = c1(k);

    if abs(A) > polynomialTolerance
        discriminant = B^2 - 4 * A * C;

        discriminantTolerance = ...
            options.polynomialTolerance ...
            * max([B^2, abs(4 * A * C), 1]);

        if discriminant >= -discriminantTolerance
            discriminant = max(discriminant, 0);
            rootDiscriminant = sqrt(discriminant);

            % Stable quadratic formula.
            if B >= 0
                q = -0.5 * (B + rootDiscriminant);
            else
                q = -0.5 * (B - rootDiscriminant);
            end

            if abs(q) > polynomialTolerance
                rootsCandidate = [q / A; C / q];
            else
                rootsCandidate = -B / (2 * A);
            end

            for j = 1:numel(rootsCandidate)
                sCandidate = rootsCandidate(j);

                if sCandidate > options.rootTolerance && ...
                        sCandidate < 1 - options.rootTolerance

                    if all(abs(candidates(1:numCandidates) - sCandidate) > ...
                            options.rootTolerance)
                        numCandidates = numCandidates + 1;
                        candidates(numCandidates) = sCandidate;
                    end

                elseif abs(sCandidate) <= options.rootTolerance
                    % Already represented by endpoint s = 0.

                elseif abs(sCandidate - 1) <= options.rootTolerance
                    % Already represented by endpoint s = 1.
                end
            end
        end

    elseif abs(c2(k)) > polynomialTolerance
        % Cubic coefficient is negligible; derivative is linear.
        sCandidate = -c1(k) / (2 * c2(k));

        if sCandidate > options.rootTolerance && ...
                sCandidate < 1 - options.rootTolerance
            numCandidates = numCandidates + 1;
            candidates(numCandidates) = sCandidate;
        end
    end

    candidates = candidates(1:numCandidates);

    candidateValues = ...
        ((c3(k) * candidates + c2(k)) .* candidates ...
        + c1(k)) .* candidates + c0(k);

    [maximumValue(k), indexMaximum] = max(candidateValues);

    sMaximum = candidates(indexMaximum);
    maximumX(k) = xLeft + h * sMaximum;

    numInteriorCandidates(k) = ...
        nnz(candidates > 0 & candidates < 1);

    if options.clipNegativeMaximum && maximumValue(k) < 0
        maximumValue(k) = 0;
        numNegativeClipped = numNegativeClipped + 1;
    end
end

info = struct();
info.numInteriorCandidates = numInteriorCandidates;
info.totalInteriorCandidates = sum(numInteriorCandidates);
info.numNegativeClipped = numNegativeClipped;

end
