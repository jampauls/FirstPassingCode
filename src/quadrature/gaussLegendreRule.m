function [nodes, weights] = gaussLegendreRule(order)
%GAUSSLEGENDRERULE Gauss--Legendre quadrature on [-1,1].
%
% MATLAB version: R2020b

validateattributes(order, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'order');

if order == 1
    nodes = 0;
    weights = 2;
    return;
end

k = (1:(order - 1)).';

beta = k ./ sqrt(4 * k.^2 - 1);

J = diag(beta, 1) ...
    + diag(beta, -1);

[V, D] = eig(J);

nodes = diag(D);

[nodes, orderIndex] = sort(nodes);
V = V(:, orderIndex);

weights = 2 * V(1, :).^2;
weights = weights(:);

end