function reconstructed = reconstructMatrixFamilyFromFeatures( ...
    report, featureDimension)
%RECONSTRUCTMATRIXFAMILYFROMFEATURES Reconstruct a matrix trajectory.
%
% Requires featureDimension basis matrices to have been stored in report.
%
% MATLAB version: R2020b

validateattributes(featureDimension, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'featureDimension');

if featureDimension > report.numBasisMatricesStored
    error('reconstructMatrixFamilyFromFeatures:BasisUnavailable', ...
        ['Requested feature dimension %d, but only %d basis matrices ', ...
         'were stored. Re-run the analysis with a larger ', ...
         'numBasisMatricesToStore.'], ...
        featureDimension, report.numBasisMatricesStored);
end

numMatrices = report.numMatrices;
r = report.matrixDimension;

reconstructed = cell(numMatrices, 1);

for k = 1:numMatrices
    Mk = zeros(r, r);

    for ell = 1:featureDimension
        Mk = Mk ...
            + report.coefficientTrajectories(ell, k) ...
            * report.basisMatrices{ell};
    end

    reconstructed{k} = 0.5 * (Mk + Mk.');
end

end