close all; clc;
%% Parameter Identification for Full-Scale Model

time_series_yaw = load("Control\Full_Scale\Control_Full_Scale_ParameterIdentification.mat");
time_series_yaw = time_series_yaw.data;

%% Plot
figure;
hold on;
grid on;

title('Yaw Vs. Time (Full-Scale Model)')
legend('show')
xlabel('Time in Seconds')
ylabel('Yaw in °')

plot(time_series_yaw, 'DisplayName','Yaw in °')


%% Finding critical gain b0

% Construct line
t1 = time_series_yaw.Time(3000);
t2 = time_series_yaw.Time(9000);
delta_t = t2 - t1;

y1 = time_series_yaw.Data(3000);
y2 = time_series_yaw.Data(9000);
delta_y = y2 - y1; 
slope =  delta_y / delta_t;
intercept = y1 - slope * t1;

t = linspace(0,60, 100);
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
u_steady_state = 750;
b_0 = delta_y/(delta_t*T*u_steady_state);
disp(['The Citrical Gain Estimate is: ', num2str(b_0)])


%% Find Dampening Coefficient a0

% For t >> T, the function approaches y(t>>T) = Ki * (t-T)*u_steady_state
% b0 is Ki/T And the dampening coefficient is 1/T 

a_1 = 1/T;
disp(['The Drag Coefficient Estimate is: ', num2str(a_1)])
%% Check by Plotting Step Response

s = tf("s");

simplified_plant = u_steady_state*b_0/((s^2+a_1*s));
figure;
step(simplified_plant,60,'r-')

hold on

plot(time_series_yaw,"DisplayName",'High Fidelity Simulated')
grid on;

title('Simplified Model vs. High Fidelity Simulation')
xlabel('Time (s)')
ylabel('Yaw (°)')
legend('show')