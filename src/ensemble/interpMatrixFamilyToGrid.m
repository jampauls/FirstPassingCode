function AInterp = interpMatrixFamilyToGrid(xSource, A, xTarget)
%INTERPMATRIXFAMILYTOGRID Interpolate a matrix family onto a new grid.
%
%   AInterp = interpMatrixFamilyToGrid(xSource, A, xTarget)
%
% Every entry of A{n} is linearly interpolated in x. xTarget must lie
% within [xSource(1), xSource(end)]; callers are responsible for clipping
% to the overlapping streamwise range beforehand.
%
% MATLAB version: R2020b

xSource = xSource(:);
xTarget = xTarget(:);

Nx = numel(xSource);

if numel(A) ~= Nx
    error('interpMatrixFamilyToGrid:LengthMismatch', ...
        'numel(A) must equal numel(xSource).');
end

if min(xTarget) < xSource(1) - 1e-9 || max(xTarget) > xSource(end) + 1e-9
    error('interpMatrixFamilyToGrid:TargetOutsideRange', ...
        'xTarget extends outside the available source range.');
end

r = size(A{1}, 1);

% Stack into a 3-D array for a single vectorized interpolation call.
Astack = zeros(r, r, Nx);
for n = 1:Nx
    Astack(:, :, n) = A{n};
end

AstackInterp = interp1( ...
    xSource, ...
    reshape(Astack, r * r, Nx).', ...
    xTarget, ...
    'linear');

AstackInterp = reshape(AstackInterp.', r, r, numel(xTarget));

Nt = numel(xTarget);
AInterp = cell(Nt, 1);

for n = 1:Nt
    An = AstackInterp(:, :, n);
    AInterp{n} = 0.5 * (An + An.');
end

end
