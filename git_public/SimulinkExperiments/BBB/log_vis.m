data_recieve = load('D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\udp_recieve_log.mat');
data_send = load('D:\AAStudium\Study\Semester_6\BScThesis\SimulinkExperiments\BBB\udp_send_log.mat');


time_send = data_send.ans(1, :);    % First row is time
value_send = data_send.ans(2, :);   % Second row is value

time_recieve = data_recieve.ans(1, :);
value_recieve = data_recieve.ans(2, :);

% Create the plot
figure;
stairs(time_send, value_send, 'Color', 'red');    
hold on;  % Ensure the next plot is on the same figure
stairs(time_recieve, value_recieve, 'Color', 'blue');
title('Value vs. Time');
xlabel('Time');
ylabel('Value');
grid on;
hold off;  % Release the plot hold
legend({'Send Data', 'Receive Data'});  % Add a legend to distinguish the lines
