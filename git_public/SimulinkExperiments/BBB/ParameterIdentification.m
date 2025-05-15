clear; close all; clc;
%% Finding b_0


measurements_angles = load("BBB\Logs\08_critgain\HENRIK_LOG_CS_Inbuilt.mat").ans;
% Reset it/ Make it start at 0° 
measurements_angles = [measurements_angles(1:3,:); (measurements_angles(4,:)-measurements_angles(4,1))];
motor_signals = load("BBB\Logs\08_critgain\HENRIK_LOG_CS_Inbuilt.mat").ans;

plot(measurements_angles(1,:),measurements_angles(4,:))
grid on;
hold on;
title('Yaw vs. Time')
xlabel('Time (s)')
ylabel('Yaw (°)')
xlim([0 measurements_angles(1,end)])

% Construct line
t1 = measurements_angles(1,end-500);
t2 = measurements_angles(1,end);
delta_t = t2 - t1;

y1 = measurements_angles(4,end-500);
y2 = measurements_angles(4,end);
delta_y = y2 - y1; 
slope =  delta_y / delta_t;
intercept = y1 - slope * t1;

t = linspace(0,t2, 10);
y = slope * t + intercept;
% Plot the line
plot(t, y, 'r--', 'LineWidth', 0.5); 
% Plot points
plot([t1, t2], [y1, y2], 'ro', 'MarkerSize', 8);

% Find T
T = -intercept/slope;
plot(T,0,'go', 'MarkerSize', 8)

plot([0 T], [0 0],'g--', 'LineWidth', 0.5, 'DisplayName','T')

% Compute Critical Gain
u_steady_state = 100;
b_0 = delta_y/(delta_t*T*u_steady_state);
disp(['The Citrical Gain Estimate is: ', num2str(b_0)])

%% Find Dampening Coefficient a0

% For t >> T, the function approaches y(t>>T) = Ki * (t-T)*u_steady_state
% b0 is Ki/T And the dampening coefficient is 1/T 

a_1 = 1/T;
disp(['The Drag Coefficient Estimate is: ', num2str(a_1)])
%% Check by Plotting Step Response

s = tf("s");

simplified_plant = 100*b_0/((s^2+a_1*s));
figure;
step(simplified_plant,35.5,'r-')

Y = measurements_angles(4,:)';
U = 100*eye(1,3545)';
hold on


plot(measurements_angles(1,:),measurements_angles(4,:),"DisplayName",'Measured')
grid on;

title('Simplified Model vs. Measurements')
xlabel('Time (s)')
ylabel('Yaw (°)')
legend('show')
xlim([0 35])

%% Get simplified state-space model

A = [0 1;
     0 -a_1];
B = [0;
    100*b_0];

C = [1 0];

D = 0;

sys_simpl = ss(A, B, C, D);

step(sys_simpl,35.5)