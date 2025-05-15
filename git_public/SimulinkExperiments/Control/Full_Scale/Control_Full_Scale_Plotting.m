%% Full-Scale Plotting

close all; clc;

% Yaw
yaw_in_deg_LQRADRC_WITH_MODEL = load("Control\Full_Scale\yaw_in_deg_LQRADRCWITHMODEL2.mat");
yaw_in_deg_LQRADRC_WITH_MODEL = yaw_in_deg_LQRADRC_WITH_MODEL.data;
yaw_in_deg_LQRADRC_NO_MODEL = load("Control\Full_Scale\yaw_in_deg_LQRADRCNOMODEL2.mat");
yaw_in_deg_LQRADRC_NO_MODEL = yaw_in_deg_LQRADRC_NO_MODEL.data;
yaw_in_deg_PLAINLQR = load("Control\Full_Scale\yaw_in_deg_LQR2.mat");
yaw_in_deg_PLAINLQR = yaw_in_deg_PLAINLQR.data;

RMSE_LQRADRC_WITH_MODEL = load("Control\Full_Scale\RMSE_LQRADRCWITHMODEL2.mat");
RMSE_LQRADRC_WITH_MODEL = RMSE_LQRADRC_WITH_MODEL.data;
RMSE_LQRADRC_NO_MODEL = load("Control\Full_Scale\RMSE_LQRADRCNOMODEL2.mat");
RMSE_LQRADRC_NO_MODEL = RMSE_LQRADRC_NO_MODEL.data;
RMSE_PLAINLQR = load("Control\Full_Scale\RMSE_LQR2.mat");
RMSE_PLAINLQR = RMSE_PLAINLQR.data;


%% Plot Yaw
figure;
hold on;
grid on;

title('Yaw Vs. Time (Full-Scale Model)')
legend('show')
xlabel('Time in Seconds')
ylabel('Yaw in °')

plot(yaw_in_deg_LQRADRC_NO_MODEL, 'DisplayName','LQR-ADRC without Model',"Color",'b',"LineWidth", 1.5)
plot(yaw_in_deg_LQRADRC_WITH_MODEL, 'DisplayName','LQR-ADRC with Model',"Color",'y',"LineWidth", 1)
plot(yaw_in_deg_PLAINLQR, 'DisplayName','Plain LQR',"Color",'r',"LineWidth", 1)

%% Plot RMSE
figure;
hold on;
grid on;

title('RSME Vs. Time (Full-Scale Model)')
legend('show')
xlabel('Time in Seconds')
ylabel('RMSE in °')

plot(RMSE_LQRADRC_NO_MODEL, 'DisplayName','LQR-ADRC without Model',"Color",'b',"LineWidth", 1.5)
plot(RMSE_LQRADRC_WITH_MODEL, 'DisplayName','LQR-ADRC with Model',"Color",'y',"LineWidth", 1)
plot(RMSE_PLAINLQR, 'DisplayName','Plain LQR',"Color",'r',"LineWidth", 1)