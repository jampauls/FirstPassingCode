function result = computeQuadraticFacetDistanceMatrix( ...
    basisMatrices, coefficientMatrix)
%COMPUTEQUADRATICFACETDISTANCEMATRIX Deterministic quadratic-form distances.
%
% For
%
%   A_k = sum_l theta_lk F_l,
%   Y_k = w.'*A_k*w,
%   w ~ N(0,I),
%
% this function computes
%
%   d_ij^2 = E[(Y_i - Y_j)^2]
%          = 2*trace((A_i-A_j)^2)
%            + trace(A_i-A_j)^2.
%
% The calculation is performed in feature coordinates:
%
%   M_lq  = trace(F_l F_q),
%   tau_l = trace(F_l),
%
%   d_ij^2
%       = 2*dtheta.'*M*dtheta
%         + (tau.'*dtheta)^2.
%
% No stochastic directions are sampled.
%
% Inputs
% ------
% basisMatrices:
%   Cell array containing symmetric feature matrices F_l.
%
% coefficientMatrix:
%   m-by-N coefficient matrix. Column k contains theta_k.
%
% Output
% ------
% result.distanceMatrix:
%   N-by-N matrix of root-mean-square quadratic-form distances.
%
% MATLAB version: R2020b

m = numel(basisMatrices);

if size(coefficientMatrix, 1) ~= m
    error('computeQuadraticFacetDistanceMatrix:DimensionMismatch', ...
        ['The coefficient matrix must have one row for each ', ...
         'feature matrix.']);
end

if m == 0
    error('computeQuadraticFacetDistanceMatrix:EmptyBasis', ...
        'At least one feature matrix is required.');
end

r = size(basisMatrices{1}, 1);

featureGram = zeros(m, m);
traceVector = zeros(m, 1);

for ell = 1:m
    Fell = basisMatrices{ell};

    if ~isequal(size(Fell), [r, r])
        error('computeQuadraticFacetDistanceMatrix:MatrixSizeMismatch', ...
            'All feature matrices must have the same square size.');
    end

    Fell = 0.5 * (Fell + Fell.');

    traceVector(ell) = trace(Fell);

    for q = 1:ell
        Fq = basisMatrices{q};
        Fq = 0.5 * (Fq + Fq.');

        value = sum(sum(Fell .* Fq));

        featureGram(ell, q) = value;
        featureGram(q, ell) = value;
    end
end

featureGram = 0.5 * (featureGram + featureGram.');

% The second-moment metric is
%
%   H = 2*M + tau*tau.'.
metricMatrix = ...
    2 * featureGram ...
    + traceVector * traceVector.';

metricMatrix = ...
    0.5 * (metricMatrix + metricMatrix.');

numConstraints = size(coefficientMatrix, 2);

% Compute squared distances through the Gram identity:
%
%   ||theta_i-theta_j||_H^2
%       = h_i + h_j - 2*theta_i.'*H*theta_j.
metricProducts = ...
    coefficientMatrix.' ...
    * metricMatrix ...
    * coefficientMatrix;

metricProducts = ...
    0.5 * (metricProducts + metricProducts.');

metricNormSquared = diag(metricProducts);

distanceSquared = ...
    bsxfun(@plus, metricNormSquared, metricNormSquared.') ...
    - 2 * metricProducts;

% Remove small negative values caused by floating-point cancellation.
distanceScale = max(max(abs(distanceSquared)), 1);

substantialNegative = ...
    distanceSquared < -1e-10 * distanceScale;

if any(substantialNegative(:))
    warning( ...
        'computeQuadraticFacetDistanceMatrix:NegativeDistanceSquared', ...
        ['The feature metric produced squared distances below the ', ...
         'expected numerical tolerance. Minimum value = %.6e.'], ...
        min(distanceSquared(:)));
end

distanceSquared = max(distanceSquared, 0);

distanceMatrix = sqrt(distanceSquared);

distanceMatrix(1:(numConstraints + 1):end) = 0;
distanceMatrix = 0.5 * (distanceMatrix + distanceMatrix.');

result = struct();

result.distanceMatrix = distanceMatrix;
result.distanceSquared = distanceSquared;

result.featureGram = featureGram;
result.traceVector = traceVector;
result.metricMatrix = metricMatrix;

result.featureDimension = m;
result.numConstraints = numConstraints;

result.isDeterministic = true;
result.interpretation = ...
    ['Root-mean-square distance between Gaussian quadratic forms: ', ...
     'sqrt(E[(Y_i-Y_j)^2]).'];

end
