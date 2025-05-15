clc; close all;
function plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval, curveNames, plotTitle)
    % PLOTCONTROLSYSTEMLOGSMULTIPLEFOLDERS Plots specific rows of data from .mat files across multiple folders,
    % allowing customization of curve names and plot title.
    % INPUT:
    %   folderPaths: Cell array of paths to folders containing .mat files.
    %   variablesToPlotPerFolder: Cell array where each folder's variables and rows are specified.
    %   plotStyles: Cell array specifying color and linewidth for each variable.
    %   timeInterval: Two-element vector specifying the desired time range, e.g., [startTime endTime].
    %   curveNames: Cell array specifying names for each plotted curve (for the legend).
    %   plotTitle: String specifying the title of the plot.

    % Initialize the plot
    figure;
    hold on;
    grid on;
    xlabel('Time');
    ylabel('Data');
    title(plotTitle); % Use the specified plot title

    styleCounter = 1; % Initialize a counter to track styles across all variables

    % Loop through each folder and its associated variables
    legendEntries = {}; % Collect legend entries
    curveCounter = 1; % Track which curve name to use from curveNames

    for folderIdx = 1:length(folderPaths)
        folderPath = folderPaths{folderIdx};
        variablesToPlot = variablesToPlotPerFolder{folderIdx};

        % Extract the folder name (last part of the folder path)
        [~, folderName] = fileparts(folderPath);

        % Loop through each variable and its associated rows
        for varIdx = 1:2:length(variablesToPlot)
            varName = variablesToPlot{varIdx};
            rowsToPlot = variablesToPlot{varIdx + 1}; % Rows to plot

            % Retrieve style for the current variable
            color = plotStyles{styleCounter}{1};
            lineWidth = plotStyles{styleCounter}{2};
            styleCounter = styleCounter + 1; % Increment the style counter

            % Construct the expected filename
            fileName = fullfile(folderPath, ['HENRIK_LOG_CS_' varName '.mat']);

            % Check if the file exists
            if isfile(fileName)
                % Load the file
                data = load(fileName);
                time = data.ans(1, :); % Extract time (1st row)

                % Apply the desired time interval
                timeIndices = (time >= timeInterval(1)) & (time <= timeInterval(2));
                filteredTime = time(timeIndices);

                % Validate rows to plot
                if max(rowsToPlot) > size(data.ans, 1)
                    warning('File %s does not contain enough rows for %s. Skipping.', fileName, varName);
                    continue;
                end

                % Extract and plot the specified rows within the time interval
                values = data.ans(rowsToPlot, timeIndices);
                % Add the curve name to the legend (if available)
                if curveCounter <= length(curveNames)
                    displayName = curveNames{curveCounter};
                else
                    displayName = strrep([folderName ' - ' varName], '_', ' '); % Fallback to folder and variable name
                end
                legendEntries{end + 1} = displayName;
                plot(filteredTime, values', 'Color', color, 'LineWidth', lineWidth, 'DisplayName', displayName); % Transpose values for proper plotting
                curveCounter = curveCounter + 1; % Move to the next curve name
            else
                warning('File %s not found.', fileName);
            end
        end
    end
    

    xline(217, '--b', 'Unmounting Flex S.', 'HandleVisibility', 'off');
    xline(219.5, '--b', '', 'HandleVisibility', 'off');
    % Add a legend to the plot
    legend('show');
    hold off;
end

% Amfitrack Roll & Pitch Vs. Flex Suppressor
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\33_FLEX_Supressor'
    
};
variablesToPlotPerFolder = {
    {'Amfitrack_Euler', [2],'Amfitrack_Euler',[3],'Amfitrack_Euler',[4]}
};
plotStyles = {
    {'g', 0.5},...
    {'b',0.5},...
    {'r',0.5}
};
timeInterval = [0, 320]; % Specify the desired time range (e.g., from 0 to 10 seconds)
curveNames = {'Amfitrack Roll Estimate', 'Amfitrack Pitch Estimate','Amfitrack Yaw Estimate'};
plotTitle = 'Amfitrack Euler Estimates during Rotation (with and without FLEX SUPPRESSOR)';

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval, curveNames, plotTitle);