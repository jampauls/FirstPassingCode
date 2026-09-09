function plotEnsembleEnergyDiagnostics(x, diagnostics, cfg)
%PLOTENSEMBLEENERGYDIAGNOSTICS Plot ensemble energy analysis outputs.

x = x(:);
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figureDirectory = fullfile(projectRoot, 'figures');
if ~isfolder(figureDirectory)
    [created, message] = mkdir(figureDirectory);
    if ~created
        error('plotEnsembleEnergyDiagnostics:FigureDirectory', ...
            'Could not create figure directory "%s": %s', ...
            figureDirectory, message);
    end
end

figure('Color', 'w', 'Name', 'Ensemble energy diagnostics');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(x, diagnostics.meanEnergyByMode, 'LineWidth', 0.75);
hold on;
modeIndex = diagnostics.mostUnstableMode;
plot(x, diagnostics.meanEnergyByMode(:, modeIndex), 'k', 'LineWidth', 2);
xlabel('x'); ylabel('Mean energy'); title('Mean energy trajectories');
grid on; box on;

nexttile;
scatter(diagnostics.beta, diagnostics.omega, 45, diagnostics.peakMeanEnergy, 'filled');
xlabel('\beta'); ylabel('\omega'); title('Peak mean energy');
colorbar; grid on; box on;

nexttile;
scatter(diagnostics.beta, diagnostics.omega, 45, diagnostics.peakGrowth, 'filled');
xlabel('\beta'); ylabel('\omega'); title('Peak growth / inlet energy');
colorbar; grid on; box on;

nexttile;
bar(diagnostics.pdf.binCenters, diagnostics.pdf.density, 1, ...
    'FaceColor', [0.2, 0.45, 0.7], 'EdgeColor', 'none');
xlabel('Energy'); ylabel('PDF estimate');
title(sprintf('PDF at x = %.6g, mode %d', ...
    diagnostics.pdf.station, diagnostics.pdf.modeIndex));
grid on; box on;

sgtitle('Ensemble energy diagnostics');

savefig(gcf, fullfile(figureDirectory, 'ensemble_energy_diagnostics.fig'));
print(gcf, fullfile(figureDirectory, 'ensemble_energy_diagnostics.png'), ...
    '-dpng', '-r150');

if nargin >= 3 && isfield(cfg, 'plot') && ...
        isfield(cfg.plot, 'closeAfterSave') && cfg.plot.closeAfterSave
    close(gcf);
end

end
