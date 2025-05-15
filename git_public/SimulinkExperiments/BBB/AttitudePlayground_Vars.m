%% Angular Measurement Playground

% Alpha for complementary filter
function alpha = calc_alpha(T,delta_t)
    alpha = (T/delta_t)/(1+T/delta_t);
end

% Standard diviation estimation for qekf
function sigma = calc_sigma(pp)
    sigma = pp/6;
end

% Kalman

delta_t = 0.01;

R = 1.8 * diag([0.001396263^2, 0.0016580628^2, 0.05^2, ((pi/180) * ...
            (0.05))^2, ((pi/180) * (0.07))^2, ((pi/180) * (0.07))^2]);
Q = diag([10^-7, 10^-7, 10^-6, 10^-8, 10^-8, 10^-7]);

A = [eye(3), delta_t * eye(3);
     zeros(3), eye(3)];

B = zeros(6,3);
C = eye(6);

% Complementary Filter

tau_roll_comp = 0.5;
tau_pitch_comp = 0.5;
tau_yaw_comp = 3;

alpha_comp_roll = calc_alpha(0.1, delta_t);
alpha_comp_pitch = calc_alpha(0.1, delta_t);
alpha_comp_yaw = calc_alpha(1, delta_t);


initial_yaw = deg2rad(0);
% Mahony
Ts = delta_t;
Ts2 = delta_t/2;
q_init = [cos(initial_yaw/2), 0, 0, sin(initial_yaw/2)]; 
k_p_mahony = 0.42; %0.42;
k_i_mahony = 0.0057; %0.0057;
Ts_IMU = delta_t;

% Bias Estimation using discrete ESO
% 1st order
N = 1;
k_ESO = 1; 
T_settle_98 = 30; % Seconds
omega_cl = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;
b0_ESO = 1;
% z-domain pole locations
z_cl = exp(-omega_cl * delta_t);
z_ESO = exp(-k_ESO * omega_cl * delta_t);
% Controller gains
k1_cESO = (1 - z_cl)/delta_t;
% Observer Gains
l1_ESO = 1 - z_ESO^2;
l2_ESO = 1/delta_t * (1 - z_ESO)^2;

l = [l1_ESO, l2_ESO]';

A_ESO = [1-l1_ESO, delta_t-l1_ESO*delta_t;
         -l2_ESO , 1 - l2_ESO*delta_t];
b_ESO = [b0_ESO*delta_t - l1_ESO*b0_ESO*delta_t;
            -l2_ESO*b0_ESO*delta_t];

% Quaternion EKF
P_0 = eye(4);

magnetic_dip_angle_theta = deg2rad(70);
r_qekf = 1/sqrt(cos(magnetic_dip_angle_theta)^2 + ...
    sin(magnetic_dip_angle_theta)^2) ...
* transpose([cos(magnetic_dip_angle_theta),0, sin(magnetic_dip_angle_theta)]);

mag_pp = 6.6; % µT
acc_pp = 2.409e-01; % m/s^2
mag_sigma = 0.8 * calc_sigma(mag_pp);
acc_sigma = 0.2 * calc_sigma(acc_pp);

% Quaternion EKF + Gyro Bias
gyro_sigma =  1*((pi/180) * (0.07));
gyro_sigma_3d = [1*gyro_sigma, 1*gyro_sigma, gyro_sigma];
bias_sigma = 0.000000001;

P_0_gb = eye(7);
bias_0 = [0, 0, 0];

% Quaternion EKF + Amfitrack
q_sigma = 0.001;

% P3022-V1- CW360 + V-Divider Circuit Parameters

voltage_divider_Vout_max = 1.728; %V

% Zero Parameters
step_time = 10;


%% %% NEW AREA %% %%

%% Lightweight Kalman

A_kalman = [1 delta_t;
            0    1  ]; 
C_kalman = eye(2);

Q_kalman = diag([0.00001, (0.07)^2]);
R_kalman = diag([(calc_sigma(0.8))^2, (0.07)^2]);
