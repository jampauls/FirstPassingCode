function plotSyntheticProblem(x, meta)
%PLOTSYNTHETICPROBLEM Plot the eigenvalue and rotation histories.

if ~isfield(meta, 'lambdaHistory') || ...
        ~isfield(meta, 'thetaHistory')
    error('plotSyntheticProblem:MissingMetadata', ...
        'The supplied metadata do not describe the synthetic problem.');
end

figure('Color', 'w', ...
    'Name', 'Synthetic transient-growth problem');

tiledlayout(2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

plot(x, meta.lambdaHistory(:, 1), ...
    'LineWidth', 2, ...
    'DisplayName', '\lambda_1(x)');

hold on;

plot(x, meta.lambdaHistory(:, 2), ...
    'LineWidth', 2, ...
    'DisplayName', '\lambda_2(x)');

xlabel('x');
ylabel('Eigenvalue of G(x)');
title('Transient normalized-energy gains');
legend('Location', 'best');
grid on;
box on;

nexttile;

plot(x, meta.thetaHistory, ...
    'LineWidth', 2);

xlabel('x');
ylabel('\theta(x)');
title('Rotation of energetic directions');
grid on;
box on;

end