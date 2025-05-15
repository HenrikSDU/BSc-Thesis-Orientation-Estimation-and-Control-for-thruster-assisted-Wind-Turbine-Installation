%% Logger & Visualization/Plotting
function plot_euler(data)
time = data(1, :);  % First row is time
curves = data(2:end, :);  % Remaining rows are the curves

% Create labels for the 24 curves
curveNames = arrayfun(@(x) sprintf('Curve %d', x), 1:24, 'UniformOutput', false);

% Plot the curves
figure;
hold on; % To overlay all plots
for i = 1:24
    plot(time, curves(i, :), 'DisplayName', curveNames{i}); % Plot each curve
end
hold off;

% Add legend and labels
legend('show');
xlabel('Time');
ylabel('Value');
title('24 Curves Plotted Against Time');

end

plot_euler(EULER_var)

