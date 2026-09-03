function result = evaluateQuadraticFeatureCgf( ...
    basisMatrices, t, computeHessian)
%EVALUATEQUADRATICFEATURECGF Evaluate the exact feature CGF.
%
% For:
%
%   X_l = w.'*F_l*w,
%   w ~ N(0,I),
%
% the cumulant-generating function is:
%
%   K(t) = -0.5*log(det(I - 2*sum_l t_l F_l)).
%
% The real CGF exists only when:
%
%   I - 2*sum_l t_l F_l
%
% is positive definite.
%
% Gradient:
%
%   dK/dt_l = trace(R*F_l),
%
% Hessian:
%
%   d2K/(dt_i dt_j)
%       = 2*trace(R*F_i*R*F_j),
%
% where:
%
%   R = (I - 2*sum_l t_l F_l)^(-1).
%
% MATLAB version: R2020b

if nargin < 3
    computeHessian = true;
end

m = numel(basisMatrices);

if m == 0
    error('evaluateQuadraticFeatureCgf:EmptyBasis', ...
        'At least one feature matrix is required.');
end

r = size(basisMatrices{1}, 1);

t = real(t(:));

if numel(t) ~= m
    error('evaluateQuadraticFeatureCgf:DimensionMismatch', ...
        't has %d entries; expected %d.', numel(t), m);
end

combinedMatrix = zeros(r, r);

for ell = 1:m
    Fell = basisMatrices{ell};

    if any(size(Fell) ~= [r, r])
        error('evaluateQuadraticFeatureCgf:FeatureSizeMismatch', ...
            'Feature matrix %d has inconsistent dimensions.', ell);
    end

    Fell = 0.5 * (Fell + Fell.');

    combinedMatrix = ...
        combinedMatrix + t(ell) * Fell;
end

combinedMatrix = ...
    0.5 * (combinedMatrix + combinedMatrix.');

domainMatrix = eye(r) - 2 * combinedMatrix;
domainMatrix = 0.5 * (domainMatrix + domainMatrix.');

[L, cholFlag] = chol(domainMatrix, 'lower');

result = struct();

result.isInDomain = cholFlag == 0;
result.t = t;

if cholFlag ~= 0
    result.K = Inf;
    result.gradient = nan(m, 1);
    result.hessian = nan(m, m);
    result.minimumDomainEigenvalue = ...
        min(real(eig(domainMatrix)));
    return;
end

logDeterminant = ...
    2 * sum(log(diag(L)));

K = -0.5 * logDeterminant;

% R = domainMatrix^(-1), evaluated through the Cholesky factors.
R = L' \ (L \ eye(r));
R = 0.5 * (R + R.');

gradient = zeros(m, 1);

RF = cell(m, 1);

for ell = 1:m
    RF{ell} = R * basisMatrices{ell};
    gradient(ell) = real(trace(RF{ell}));
end

if computeHessian
    hessian = zeros(m, m);

    for i = 1:m
        for j = i:m
            value = 2 * real(trace( ...
                RF{i} * RF{j}));

            hessian(i, j) = value;
            hessian(j, i) = value;
        end
    end

    hessian = 0.5 * (hessian + hessian.');
else
    hessian = [];
end

result.K = K;
result.gradient = gradient;
result.hessian = hessian;

result.domainMatrix = domainMatrix;
result.inverseDomainMatrix = R;

result.minimumDomainEigenvalue = ...
    min(real(eig(domainMatrix)));

end