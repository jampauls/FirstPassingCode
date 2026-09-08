function plotTransitionResults(x, prob, pLocal, post, cfg)
%PLOTTRANSITIONRESULTS Plot principal first-transition outputs.

x = x(:);

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figureDirectory = fullfile(projectRoot, 'figures');

if ~isfolder(figureDirectory)
    [created, message] = mkdir(figureDirectory);
    if ~created
        error('plotTransitionResults:FigureDirectory', ...
            'Could not create figure directory "%s": %s', ...
            figureDirectory, message);
    end
end

figure('Color', 'w', ...
    'Name', 'First-transition probability');

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% CDF and local exceedance
% -------------------------------------------------------------------------
nexttile;
hold on;

plot(x, prob.F, ...
    'LineWidth', 2, ...
    'DisplayName', 'First-transition CDF');

if ~isempty(pLocal)
    plot(x, pLocal, '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Local exceedance');
end

if isfield(prob, 'standardError') && ...
        any(prob.standardError > 0)

    lowerBand = max( ...
        prob.F - 2 * prob.standardError, 0);

    upperBand = min( ...
        prob.F + 2 * prob.standardError, 1);

    fill( ...
        [x; flipud(x)], ...
        [lowerBand; flipud(upperBand)], ...
        [0.75, 0.85, 1.00], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, ...
        'DisplayName', 'Approx. 2-SE RQMC band');

    % Redraw the mean over the shaded region.
    plot(x, prob.F, ...
        'LineWidth', 2, ...
        'DisplayName', 'First-transition CDF');
end


xlabel('x');
ylabel('Probability');
title('Transition probability');
grid on;
box on;
ylim([0, 1]);
legend('Location', 'best');

% -------------------------------------------------------------------------
% Survival probability
% -------------------------------------------------------------------------
nexttile;

plot(x, prob.S, ...
    'LineWidth', 2);

xlabel('x');
ylabel('S(x)');
title('Survival probability');
grid on;
box on;
ylim([0, 1]);

% -------------------------------------------------------------------------
% Interval probabilities
% -------------------------------------------------------------------------
nexttile;

if isfield(cfg.plot, 'showIntervalProbabilities') && ...
        cfg.plot.showIntervalProbabilities

    stem(x, prob.pInterval, ...
        'filled', ...
        'MarkerSize', 3);

    xlabel('x');
    ylabel('Interval probability');
    title('First-transition mass by interval');
    grid on;
    box on;
else
    axis off;
end

% -------------------------------------------------------------------------
% Memory correction
% -------------------------------------------------------------------------
nexttile;

if isfield(post, 'memoryCorrection') && ...
        ~isempty(post.memoryCorrection)

    hold on;

    if isfield(post, 'memoryCorrectionStandardError') && ...
            ~isempty(post.memoryCorrectionStandardError)

        correctionLower = ...
            max( ...
            post.memoryCorrection ...
            - 2 * post.memoryCorrectionStandardError, ...
            0);

        correctionUpper = ...
            post.memoryCorrection ...
            + 2 * post.memoryCorrectionStandardError;

        fill( ...
            [x; flipud(x)], ...
            [correctionLower; flipud(correctionUpper)], ...
            [0.85, 0.85, 0.85], ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, ...
            'DisplayName', 'Approx. 2-SE band');
    end

    plot(x, post.memoryCorrection, ...
        'LineWidth', 2, ...
        'DisplayName', 'Memory correction');

    xlabel('x');
    ylabel('F_{X_{tr}}(x) - p_{local}(x)');
    title('Memory correction');
    grid on;
    box on;

else
    axis off;
end

sgtitle('Memory-aware first-transition prediction');

% Print quantiles.
if isfield(post, 'quantiles')
    fprintf('\nTransition quantiles:\n');

    for j = 1:numel(post.quantiles.levels)
        level = post.quantiles.levels(j);

        if post.quantiles.isCensored(j)
            fprintf('  p = %.2f: not reached by x_max\n', level);
        else
            fprintf('  p = %.2f: x = %.6g\n', ...
                level, post.quantiles.values(j));
        end
    end

    fprintf('\n');
end

figureFile = fullfile(figureDirectory, 'first_transition_results.fig');
pngFile = fullfile(figureDirectory, 'first_transition_results.png');

savefig(gcf, figureFile);
print(gcf, pngFile, '-dpng', '-r150');

fprintf('Saved figure: %s\n', figureFile);
fprintf('Saved figure: %s\n', pngFile);

end