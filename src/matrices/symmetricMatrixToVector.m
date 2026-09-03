function vector = symmetricMatrixToVector(matrix)
%SYMMETRICMATRIXTOVECTOR Isometric vectorization of a symmetric matrix.
%
% For a real symmetric r-by-r matrix A, return its upper-triangular
% entries with off-diagonal entries multiplied by sqrt(2). Then:
%
%   svec(A)'*svec(B) = trace(A*B)
%
% for real symmetric A and B.
%
% MATLAB version: R2020b

if ~isnumeric(matrix) || ndims(matrix) ~= 2
    error('symmetricMatrixToVector:InvalidInput', ...
        'Input must be a numeric matrix.');
end

[r1, r2] = size(matrix);

if r1 ~= r2
    error('symmetricMatrixToVector:NonSquareMatrix', ...
        'Input matrix must be square.');
end

if ~isreal(matrix)
    error('symmetricMatrixToVector:ComplexMatrix', ...
        'This function requires a real symmetric matrix.');
end

matrix = 0.5 * (matrix + matrix.');

upperMask = triu(true(r1));

[rowIndex, columnIndex] = find(upperMask);

linearIndex = sub2ind( ...
    [r1, r1], rowIndex, columnIndex);

vector = matrix(linearIndex);

offDiagonal = rowIndex ~= columnIndex;

vector(offDiagonal) = ...
    sqrt(2) * vector(offDiagonal);

vector = vector(:);

end
