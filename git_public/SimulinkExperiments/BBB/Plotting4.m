clear; close all; clc;

%% Loading Data
% Specify the folder containing the .mat files
folderPath = 'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\05_fast_varmot'; % 5_fast_varmot

% Get a list of all .mat files in the folder
fileList = dir(fullfile(folderPath, '*.mat'));

% Loop through each file and load it
for i = 1:length(fileList)
    fileName = fullfile(folderPath, fileList(i).name);
    loadedData = load(fileName);
    disp(['Loaded: ', fileList(i).name]);
    
    % Dynamically assign variables (use with caution)
    fieldNames = fieldnames(loadedData);
    for j = 1:numel(fieldNames)
        assignin('base', fieldNames{j}, loadedData.(fieldNames{j}));
    end
end

%% Accumulating and Resetting AMFITRACK_EULER
resetTime = 150; % Define the timestep for resetting accumulation
time = AMFITRACK_EULER(1, :); % Extract time from the first row
data = AMFITRACK_EULER(2:end, :); % Extract data (remaining rows)

% Initialize accumulated data and previous values
accumulativeData = zeros(size(data));
previousAccumulatedValue = zeros(size(data, 1), 1); % To store previous accumulated values

for t = 1:length(time)
    if time(t) >= resetTime
        previousAccumulatedValue = zeros(size(data, 1), 1); % Reset accumulation at the specified timestep
    end
    
    % Calculate delta for each row and adjust for wraparound (-180° to 180°)
    if t > 1 % Ensure previous timestep exists
        deltaData = data(:, t) - data(:, t - 1);
        deltaData(deltaData > 180) = deltaData(deltaData > 180) - 360;
        deltaData(deltaData < -180) = deltaData(deltaData < -180) + 360;
        previousAccumulatedValue = previousAccumulatedValue + deltaData;
    end
    
    accumulativeData(:, t) = previousAccumulatedValue;
end

% Replace original data with the accumulated data
AMFITRACK_EULER(2:end, :) = accumulativeData;
% Offset RAW [3] by 18 degrees
RAW(4, :) = RAW(4, :) + 12; % Adjusting the index to access the third row of data (MATLAB starts indexing from 1)
%% Plotting Function
%% Plotting Function
function plotData(varargin)
    % Function to plot selected rows from one or more n_rows x m_columns matrices.
    % Allows specifying rows for each signal individually. If no rows are
    % specified, all data rows (except the time row) are plotted.
    %
    % Example usage:
    % plotData(IMU_AMG, [4:7], QEKF, [3], [10, 180], [-185, 185], 'Custom Title');

    if nargin < 1
        error('At least one input is required.');
    end

    % Initialize optional inputs
    customTitle = 'Fusion: Plotting Experimental Results'; % Default title
    timeInterval = []; % Default (no restriction on time interval)
    yLimits = []; % Default (auto y-axis limits)

    % Parse inputs
    parsedSignals = {}; % Stores signal matrices
    parsedRows = {}; % Stores row selections for each signal
    idx = 1;

    while idx <= nargin
        if isnumeric(varargin{idx}) && length(varargin{idx}) == 2
            % Check if this is a y-axis limits or time interval
            if isempty(yLimits)
                yLimits = varargin{idx}; % Assign to y-axis limits
            else
                timeInterval = varargin{idx}; % Assign to time interval
            end
            idx = idx + 1;
        elseif ischar(varargin{idx}) || isstring(varargin{idx})
            customTitle = varargin{idx}; % Extract the custom title
            idx = idx + 1;
        elseif isnumeric(varargin{idx}) || ismatrix(varargin{idx})
            % This is a signal matrix
            parsedSignals{end + 1} = varargin{idx};
            if idx + 1 <= nargin && isnumeric(varargin{idx + 1})
                % The next argument is a row selection
                parsedRows{end + 1} = varargin{idx + 1};
                idx = idx + 2;
            else
                % No row selection provided, plot all rows
                parsedRows{end + 1} = [];
                idx = idx + 1;
            end
        else
            error('Invalid input format.');
        end
    end

    % Create the figure and start plotting
    figure;
    hold on;

    for k = 1:length(parsedSignals)
        inputData = parsedSignals{k};
        rowSelection = parsedRows{k};

        % Extract the variable name (for DisplayName)
        variableName = inputname(k * 2 - 1); % Adjusted for function arguments
        if isempty(variableName)
            variableName = ['Variable ', num2str(k)];
        end
        variableName = strrep(variableName, '_', ' '); % Replace underscores with spaces

        % Check if the matrix has at least two rows (time row + data rows)
        if size(inputData, 1) < 2
            error(['Matrix must have at least two rows: time row and data rows (', variableName, ')']);
        end

        % Extract time (first row) and data (remaining rows)
        time = inputData(1, :);
        data = inputData(2:end, :); % Skip the time row

        % Apply time interval restriction, if specified
        if ~isempty(timeInterval)
            timeMask = (time >= timeInterval(1)) & (time <= timeInterval(2));
            time = time(timeMask);
            data = data(:, timeMask);
        end

        % Plot the selected rows (or all rows if none specified)
        if isempty(rowSelection)
            rowSelection = 1:size(data, 1); % Default to all data rows
        end

        for i = rowSelection
            if i > size(data, 1)
                warning(['Row ', num2str(i), ' exceeds the number of available rows in ', variableName]);
                continue;
            end
            plot(time, data(i, :), 'DisplayName', [variableName, ' Row ', num2str(i)]);
        end
    end

    % Add labels, legend, title, and custom axis limits
    xlabel('Time');
    ylabel('Value');
    title(customTitle);
    legend('show');

    % Add horizontal reference lines without showing them in the legend
    yline(-180, '--b', '-180°', 'HandleVisibility', 'off');
    yline(-90, '--b', '-90°', 'HandleVisibility', 'off');
    yline(0, '--b', '0°', 'HandleVisibility', 'off');
    yline(90, '--b', '90°', 'HandleVisibility', 'off');
    yline(180, '--b', '180°', 'HandleVisibility', 'off');

    grid on;

    % Apply y-axis limits, if specified
    if ~isempty(yLimits)
        ylim(yLimits);
    end

    hold off;
end
plotData(QEKF_AOY, [1:3],AMFITRACK_EULER,3,MOTOR_SIGNAL,1,[-185, 185], [10, 160], 'Inbuilt vs. QEKF (Motor On)');
%plotData(IMU_AMG, [4:6], RAW, [3], INBUILT,[1:3],COMPLENTARY,3,[-180, 180], [20, 185], 'Drift in Inbuilt Sensor Fusion Motor Off');
