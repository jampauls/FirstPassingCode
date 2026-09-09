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
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(x, diagnostics.meanEnergyByMode, 'LineWidth', 0.75);
hold on;
modeIndex = diagnostics.mostUnstableMode;
plot(x, diagnostics.meanEnergyByMode(:, modeIndex), 'k', 'LineWidth', 2);
xlabel('x'); ylabel('Mean energy'); title('Mean energy trajectories');
grid on; box on;

nexttile;
plotEnergyMap(diagnostics.beta, diagnostics.omega, ...
    diagnostics.peakMeanEnergy, '\beta', '\omega', ...
    'Peak mean energy');

nexttile;
plotEnergyMap(diagnostics.beta, diagnostics.omega, ...
    diagnostics.peakGrowth, '\beta', '\omega', ...
    'Peak growth / inlet energy');

nexttile;
plotEnergyMap(x, diagnostics.omega, diagnostics.meanEnergyByMode, ...
    'x', '\omega', 'Mean energy versus x and \omega');

nexttile;
plotEnergyMap(x, diagnostics.beta, diagnostics.meanEnergyByMode, ...
    'x', '\beta', 'Mean energy versus x and \beta');

nexttile;
plotEnergyMap(diagnostics.beta, diagnostics.omega, ...
    diagnostics.inletEnergy, '\beta', '\omega', ...
    'Initial mean energy');

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

    function plotEnergyMap(horizontalValues, verticalValues, values, ...
            horizontalLabel, verticalLabel, plotTitle)

        horizontalValues = horizontalValues(:);
        verticalValues = verticalValues(:);
        values = double(values);

        if ~isvector(values) && ...
                size(values, 1) == numel(horizontalValues) && ...
                size(values, 2) == numel(verticalValues)
            values = values.';
        end

        if ~isvector(values) && ...
                size(values, 1) == numel(verticalValues) && ...
                size(values, 2) == numel(horizontalValues)
            [uniqueHorizontal, ~, horizontalIndex] = unique(horizontalValues);
            [uniqueVertical, ~, verticalIndex] = unique(verticalValues);
            [horizontalIndexGrid, verticalIndexGrid] = meshgrid( ...
                horizontalIndex, verticalIndex);

            gridValues = accumarray( ...
                [verticalIndexGrid(:), horizontalIndexGrid(:)], ...
                values(:), ...
                [numel(uniqueVertical), numel(uniqueHorizontal)], ...
                @mean, NaN);

            [horizontalGrid, verticalGrid] = meshgrid( ...
                uniqueHorizontal, uniqueVertical);
        else
            uniqueHorizontal = unique(horizontalValues);
            uniqueVertical = unique(verticalValues);
            [horizontalGrid, verticalGrid] = meshgrid( ...
                uniqueHorizontal, uniqueVertical);

            gridValues = griddata(horizontalValues, verticalValues, values(:), ...
                horizontalGrid, verticalGrid, 'nearest');
        end

        if numel(uniqueHorizontal) < 2 || numel(uniqueVertical) < 2
            plot(horizontalValues, values, 'o-');
        else
            contourf(horizontalGrid, verticalGrid, gridValues, 20, ...
                'LineColor', 'none');
        end

        xlabel(horizontalLabel);
        ylabel(verticalLabel);
        title(plotTitle);
        colorbar;
        grid on;
        box on;

    end

end
