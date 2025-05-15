close all;
function plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval)
    % PLOTCONTROLSYSTEMLOGSMULTIPLEFOLDERS Plots specific rows of data from .mat files across multiple folders,
    % with an option to specify a desired time interval.
    % INPUT:
    %   folderPaths: Cell array of paths to folders containing .mat files.
    %   variablesToPlotPerFolder: Cell array where each folder's variables and rows are specified.
    %   plotStyles: Cell array specifying color and linewidth for each variable.
    %   timeInterval: Two-element vector specifying the desired time range, e.g., [startTime endTime].

    % Initialize the plot
    figure;
    hold on;
    grid on;
    xlabel('Time');
    ylabel('Data');
    title('Control System Logs Across Different Runs');

    styleCounter = 1; % Initialize a counter to track styles across all variables

    % Loop through each folder and its associated variables
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
                % Add the folder name and variable name to the legend
                displayName = strrep([folderName ' - ' varName], '_', ' '); % Replace underscores with spaces
                plot(filteredTime, values', 'Color', color, 'LineWidth', lineWidth, 'DisplayName', displayName); % Transpose values for proper plotting
            else
                warning('File %s not found.', fileName);
            end
        end
    end

    % Add a legend to the plot
    legend('show');
    hold off;
end
%{
% Example usage:
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\REF', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\23_LADRC_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\24_LQRADRCNOMODEL_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\25_LQRADRCNOMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\31_LQRADRCWMODEL_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\30_LQRADRCWMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\28_PLAINLQR_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\29_PLAINLQR_WithDisturbance'
};
variablesToPlotPerFolder = {
    {'REF', [2]},...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}
};
plotStyles = {
    {'g', 2}, ...
    {'b', 1}, ...
    {'b', 1}, ...
    {'y', 1}, ...
    {'y', 1}, ...
    {'r', 1}, ...
    {'r', 1}, ...
    {'c', 1},...
    {'c', 1},...
    {'b', 1}, ...
    {'b', 1}, ...
    {'r', 1}, ...
    {'r', 1}, ...
    {'r', 1}, ...
    {'r', 1} 
};
timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);


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
plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);
%}
%% New runs of controllers
 
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\REF', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\34_LADRC_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\36_LQRADRCNOMODEL_NoDistubance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\37_LQRADRCNOMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\38_LQRADRCWITHMODEL_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\39_LQRADRCWITHMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\40_LQR_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\41_LQR_WithDisturbance'
};
variablesToPlotPerFolder = {
    {'REF', [2]},...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}
};
plotStyles = {
    {'g', 2}, ... % Green
    {'b', 1}, ... % Blue
    {'m', 1}, ... % Magenta (was blue)
    {'y', 1}, ... % Yellow
    {'k', 1}, ... % Black (was yellow)
    {'r', 1}, ... % Red
    {'c', 1}, ... % Cyan
    {'m', 1}, ... % Magenta
    {'y', 1}, ... % Yellow (was red)
    {'c', 1}, ... % Cyan (was red)
};

timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);

%% Better Plotting - No Diturbance

%% New runs of controllers
 
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\REF', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\34_LADRC_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\36_LQRADRCNOMODEL_NoDistubance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\38_LQRADRCWITHMODEL_NoDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\40_LQR_NoDisturbance', ...
};
variablesToPlotPerFolder = {
    {'REF', [2]},...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}
    };
plotStyles = {
    {'g', 2}, ... % Green
    {'b', 1}, ... % Blue
    {'m', 1}, ... % Magenta (was blue)
    {'y', 1}, ... % Yellow
    {'r', 1}, ... % Red
};

timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);

%% Better Plotting - With Diturbance

%% New runs of controllers
 
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\REF', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\35_LADRC_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\37_LQRADRCNOMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\39_LQRADRCWITHMODEL_WithDisturbance', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\41_LQR_WithDisturbance', ...
};
variablesToPlotPerFolder = {
    {'REF', [2]},...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}, ...
    {'QEKF_AOY', [4]}
    };
plotStyles = {
    {'g', 2}, ... % Green
    {'b', 1}, ... % Blue
    {'m', 1}, ... % Magenta (was blue)
    {'y', 1}, ... % Yellow
    {'r', 1}, ... % Black (was yellow)
};

timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);

%% Better Plotting - RMSE

%% New runs of controllers
 
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\45_LADRC_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\45_LADRC_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\42_LQRADRCNOMODEL_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\43_LQRADRCWITHMODEL_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\44_LQR_RMSE90', ...
};

variablesToPlotPerFolder = {
    {'REF', [2]},...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    {'QEKF_AOY', [4],'xhat_RMSE',[2]}, ...
    };
plotStyles = {
    {'g', 2}, ... % Green
    {'b', 1}, ... % Blue
    {'b', 1}, ... % Blue
    {'m', 1}, ... % Magenta 
    {'m', 1}, ... % Magenta 
    {'y', 1}, ... % Yellow
    {'y', 1}, ... % Yellow
    {'r', 1}, ... % Red
    {'r', 1}, ... % Red
};

timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);

%% New runs of controllers
 
folderPaths = {
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\45_LADRC_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\42_LQRADRCNOMODEL_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\43_LQRADRCWITHMODEL_RMSE90', ...
    'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\44_LQR_RMSE90', ...
};

variablesToPlotPerFolder = {
    {'xhat_RMSE',[2]}, ...
    {'xhat_RMSE',[2]}, ...
    {'xhat_RMSE',[2]}, ...
    {'xhat_RMSE',[2]}, ...
    };
plotStyles = {
    {'b', 1}, ... % Blue
    {'m', 1}, ... % Magenta 
    {'y', 1}, ... % Yellow
    {'r', 1}, ... % Red
};

timeInterval = [20, 100]; % Specify the desired time range (e.g., from 0 to 10 seconds)

plotControlSystemLogsMultipleFolders(folderPaths, variablesToPlotPerFolder, plotStyles, timeInterval);
