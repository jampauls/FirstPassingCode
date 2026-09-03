function matrix = symmetricVectorToMatrix(vector, dimension)
%SYMMETRICVECTORTOMATRIX Inverse of symmetricMatrixToVector.
%
% MATLAB version: R2020b

validateattributes(dimension, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    mfilename, 'dimension');

vector = vector(:);

expectedLength = ...
    dimension * (dimension + 1) / 2;

if numel(vector) ~= expectedLength
    error('symmetricVectorToMatrix:LengthMismatch', ...
        'Vector has length %d; expected %d.', ...
        numel(vector), expectedLength);
end

upperMask = triu(true(dimension));

[rowIndex, columnIndex] = find(upperMask);

offDiagonal = rowIndex ~= columnIndex;

values = vector;
values(offDiagonal) = ...
    values(offDiagonal) / sqrt(2);

matrix = zeros(dimension, dimension);

linearUpper = sub2ind( ...
    [dimension, dimension], ...
    rowIndex, columnIndex);

matrix(linearUpper) = values;

% Copy the strict upper triangle to the lower triangle.
strictUpper = triu(matrix, 1);

matrix = ...
    diag(diag(matrix)) ...
    + strictUpper ...
    + strictUpper.';

end
