function diagnostics = computeEnsembleEnergyDiagnostics( ...
    x, meanEnergyByMode, omega, beta, fileList, interpolatedAByMode, cfg)
%COMPUTEENSEMBLEENERGYDIAGNOSTICS Summarize per-mode energy growth.

x = x(:);
meanEnergyByMode = double(meanEnergyByMode);
M = size(meanEnergyByMode, 2);

if size(meanEnergyByMode, 1) ~= numel(x) || ...
        numel(omega) ~= M || numel(beta) ~= M || numel(fileList) ~= M
    error('computeEnsembleEnergyDiagnostics:DimensionMismatch', ...
        'Mode metadata and mean-energy trajectories are inconsistent.');
end

inletEnergy = meanEnergyByMode(1, :).';
peakMeanEnergy = max(meanEnergyByMode, [], 1).';
[~, peakStationIndexByMode] = max(meanEnergyByMode, [], 1);
peakGrowth = peakMeanEnergy ./ max(inletEnergy, eps);

[maximumPeakGrowth, mostUnstableMode] = max(peakGrowth);
mostUnstablePeakStation = peakStationIndexByMode(mostUnstableMode);

diagnostics = struct();
diagnostics.meanEnergyByMode = meanEnergyByMode;
diagnostics.inletEnergy = inletEnergy;
diagnostics.peakMeanEnergy = peakMeanEnergy;
diagnostics.peakStationIndexByMode = peakStationIndexByMode;
diagnostics.peakStationByMode = x(peakStationIndexByMode);
diagnostics.peakGrowth = peakGrowth;
diagnostics.omega = omega(:);
diagnostics.beta = beta(:);
diagnostics.fileList = fileList;
diagnostics.mostUnstableMode = mostUnstableMode;
diagnostics.mostUnstablePeakGrowth = maximumPeakGrowth;
diagnostics.mostUnstablePeakStationIndex = mostUnstablePeakStation;
diagnostics.mostUnstablePeakStation = x(mostUnstablePeakStation);

peakMatrix = interpolatedAByMode{mostUnstableMode}{mostUnstablePeakStation};
peakMatrix = 0.5 * (peakMatrix + peakMatrix.');
eigenvalues = real(eig(peakMatrix));

numSamples = 100000;
seed = 1;
if isfield(cfg, 'analysis')
    if isfield(cfg.analysis, 'pdfNumSamples')
        numSamples = cfg.analysis.pdfNumSamples;
    end
    if isfield(cfg.analysis, 'pdfRandomSeed')
        seed = cfg.analysis.pdfRandomSeed;
    end
end

rng(seed, 'twister');
sampleBlockSize = 10000;
energySamples = zeros(numSamples, 1);
for firstSample = 1:sampleBlockSize:numSamples
    lastSample = min(firstSample + sampleBlockSize - 1, numSamples);
    gaussianSamples = randn(numel(eigenvalues), lastSample - firstSample + 1);
    energySamples(firstSample:lastSample) = ...
        sum(eigenvalues .* (gaussianSamples .^ 2), 1).';
end

numBins = 100;
if isfield(cfg, 'analysis') && isfield(cfg.analysis, 'pdfNumBins')
    numBins = cfg.analysis.pdfNumBins;
end

binEdges = linspace(0, max(energySamples) * 1.05, numBins + 1).';
binWidth = binEdges(2) - binEdges(1);
binCenters = 0.5 * (binEdges(1:end-1) + binEdges(2:end));
pdfCounts = histcounts(energySamples, binEdges).';
pdfDensity = pdfCounts / (numSamples * binWidth);

diagnostics.pdf = struct();
diagnostics.pdf.energySamples = energySamples;
diagnostics.pdf.eigenvalues = eigenvalues;
diagnostics.pdf.binEdges = binEdges;
diagnostics.pdf.binCenters = binCenters;
diagnostics.pdf.density = pdfDensity;
diagnostics.pdf.numSamples = numSamples;
diagnostics.pdf.randomSeed = seed;
diagnostics.pdf.station = diagnostics.mostUnstablePeakStation;
diagnostics.pdf.modeIndex = mostUnstableMode;

end
