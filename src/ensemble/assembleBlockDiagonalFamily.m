function [Gtotal, blockSizes] = assembleBlockDiagonalFamily(GredByMode)
%ASSEMBLEBLOCKDIAGONALFAMILY Assemble a sparse block-diagonal G family.
%
%   [Gtotal, blockSizes] = assembleBlockDiagonalFamily(GredByMode)
%
% GredByMode{i} is the cell array of reduced matrices for mode i, one per
% streamwise station. Modes are assumed independent, so the total
% quadratic form is block-diagonal across the concatenated stochastic
% vector. Blocks are stored sparsely: nonzero storage scales with
% sum(k_i^2), not with (sum k_i)^2.
%
% MATLAB version: R2020b

M = numel(GredByMode);

if M == 0
    error('assembleBlockDiagonalFamily:EmptyEnsemble', ...
        'GredByMode must contain at least one mode.');
end

Nx = numel(GredByMode{1});

blockSizes = zeros(M, 1);
for i = 1:M
    if numel(GredByMode{i}) ~= Nx
        error('assembleBlockDiagonalFamily:LengthMismatch', ...
            'Every mode must contribute the same number of stations.');
    end
    blockSizes(i) = size(GredByMode{i}{1}, 1);
end

Gtotal = cell(Nx, 1);

for n = 1:Nx
    blocks = cell(M, 1);

    for i = 1:M
        blocks{i} = sparse(GredByMode{i}{n});
    end

    Gtotal{n} = blkdiag(blocks{:});
end

end
