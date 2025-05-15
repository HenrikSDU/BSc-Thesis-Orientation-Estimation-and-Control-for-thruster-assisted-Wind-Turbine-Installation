clear; close all; clc;
%% Loading Data
% Specify the folder containing the .mat files 
folderPath = 'D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\Logs\03'; % 5_fast_varmot 

% Get a list of all .mat files in the folder 
fileList = dir(fullfile(folderPath, '*.mat')); 

% Loop through each file and load it 
 for i = 1:length(fileList) 
     fileName = fullfile(folderPath, fileList(i).name); 
     loadedData = load(fileName); 
     % Load the .mat file 
     % Process or use the loaded data here 
     disp(['Loaded: ', fileList(i).name]); 
     % Example: assign variables dynamically (use with caution) 
     fieldNames = fieldnames(loadedData); 
     for j = 1:numel(fieldNames) 
     assignin('base', fieldNames{j}, loadedData.(fieldNames{j})); 
     end 
 end
%% Plotting
function plotData(varargin)

    % Function to plot data from one or more n_rows x m_columns matrices.
    % Matrices can be combined into the same plot using brackets.

    if nargin < 1
        error('At least one input is required.');
    end

    % Initialize optional inputs
    customTitle = 'Fusion: Plotting Experimental Results'; % Default title
    timeInterval = []; % Default (no restriction on time interval)
    yLimits = []; % Default (auto y-axis limits)

    % Parse inputs
    if ischar(varargin{end}) || isstring(varargin{end})
        customTitle = varargin{end}; % Extract the custom title
        varargin(end) = []; % Remove the title from the input arguments
    end

    if isnumeric(varargin{end}) && length(varargin{end}) == 2
        yLimits = varargin{end}; % Extract y-axis limits
        varargin(end) = []; % Remove from input arguments
    end

    if isnumeric(varargin{end}) && length(varargin{end}) == 2
        timeInterval = varargin{end}; % Extract time interval
        varargin(end) = []; % Remove from input arguments
    end

    figure; % Create a single figure for grouped plots
    hold on; % Allow multiple plots on the same figure
    
    for k = 1:length(varargin)
        % Get the current input, whether it's a single matrix or grouped
        inputData = varargin{k};
        
        % If input is grouped (horizontal concatenation), split into cells
        if iscell(inputData)
            inputData = cell2mat(inputData);
        end
        
        % Extract the variable name and replace underscores with spaces
        variableName = inputname(k);
        if isempty(variableName)
            variableName = ['Variable ', num2str(k)];
        end
        variableName = strrep(variableName, '_', ' '); % Replace underscores with spaces
        
        % Check if the matrix has at least two rows (time + data)
        if size(inputData, 1) < 2
            error('Matrix must have at least two rows: time row and data rows.');
        end

        % Extract time (first row) and data (remaining rows)
        time = inputData(1, :); % Time (first row)
        data = inputData(2:end, :); % Data (remaining rows)

        % Apply time interval restriction, if specified
        if ~isempty(timeInterval)
            timeMask = (time >= timeInterval(1)) & (time <= timeInterval(2));
            time = time(timeMask);
            data = data(:, timeMask);
        end

        % Plot each curve in the current matrix
        for i = 1:size(data, 1)
            plot(time, data(i, :), 'DisplayName', [variableName, ' Row ', num2str(i)]);
        end
    end
    
    % Add labels, legend, title, and custom axis limits
    xlabel('Time');
    ylabel('Value');
    title(customTitle); % Use the dynamic title
    legend('show');
    
    % Add horizontal lines without showing them in the legend
    yline(-180, '--b','-180°', 'HandleVisibility','off');
    yline(-90, '--b','-90°', 'HandleVisibility','off');
    yline(0, '--b','0°', 'HandleVisibility','off');
    yline(90, '--b','90°', 'HandleVisibility','off');
    yline(180, '--b','180°', 'HandleVisibility','off');

    grid on;

    % Apply y-axis limits, if specified
    if ~isempty(yLimits)
        ylim(yLimits);
    end

    hold off; % Stop adding to the current figure
end

plotData(INBUILT,AMFITRACK_EULER, MOTOR_SIGNAL,[10, 120], [-185,185], 'Drift in Inbuilt Sensor Fusion Motor Off');

