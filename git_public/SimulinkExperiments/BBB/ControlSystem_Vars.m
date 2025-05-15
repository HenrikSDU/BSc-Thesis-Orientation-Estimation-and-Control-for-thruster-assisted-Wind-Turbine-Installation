clear; close all; clc;
%% %%%%%%%%%%%%%%%%% Control System Variables %%%%%%%%%%%%%%%%% %%

% Sample Time
delta_t = 0.01;

% Amfitrack Settle Time
amfitrack_time = 1; % Seconds

% Perfomance Evaluation Reset (To let things settle)
reset_time = 20; % Seconds

%% Motor Mixing

motor_bias = 30;

%{
d_motors = 0.21;

MIX = [-3/(2*d_motors);
        3/(2*d_motors)];
%}
%% Sensor Fusion Parameters

% Standard diviation estimation for QEKF
function sigma = calc_sigma(pp)
    sigma = pp/6;
end

% IMU Noise
acc_pp = 2.409e-01; % m/s^2
acc_sigma = 0.2 * calc_sigma(acc_pp);
gyro_sigma =  1*((pi/180) * (0.07));
gyro_sigma_3d = [1*gyro_sigma, 1*gyro_sigma, gyro_sigma];

% Amfitrack Noise
q_sigma = 0.001; %calc_sigma(0.5);% Best working value 0.001

% Initials
P_0 = eye(4);

%% Control Parameters

% %%%%% Shared Parameters %%%%% %

% Model Parameters
a1 = 0.16072; % Drag Coefficient Estimate obtained in ParameterIdentification.m
b0 =  1.025; % Critical Gain Obtained in ParameterIdentification.m

% %%%%%%%%%%%%%%% DLADRC %%%%%%%%%%%%%%% %

N = 2;
k_ESO = 7; 
T_settle_98 = 3; % Seconds
omega_cl = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

% z-domain pole locations
z_cl = exp(-omega_cl * delta_t);
z_ESO = exp(-k_ESO * omega_cl * delta_t);

% Controller gains
k1_c_d = (1 - z_cl)^2/delta_t^2;
k2_c_d = (4-(1+z_cl)^2)/(2*delta_t);
k_T_d = [k1_c_d, k2_c_d];

% Observer Gains
l1_ESO_d = 1 - z_ESO^3;
l2_ESO_d = 3/(2*delta_t) * (1 - z_ESO)^2*(1 + z_ESO);
l3_ESO_d = 1/(delta_t^2) * (1-z_ESO)^3;
l_ESO_d = [l1_ESO_d, l2_ESO_d, l3_ESO_d]';

A_ESO_d = [1-l1_ESO_d, delta_t-l1_ESO_d*delta_t, 1/2 * delta_t^2 - 1/2 * l1_ESO_d*delta_t^2;
         -l2_ESO_d , 1 - l2_ESO_d*delta_t, delta_t - 1/2 * l2_ESO_d*delta_t^2;
        -l3_ESO_d, -l3_ESO_d*delta_t, 1-1/2*l3_ESO_d*delta_t^2;];
b_ESO_d = [1/2*b0*delta_t^2 - l1_ESO_d*b0*delta_t^2;
            b0*delta_t-1/2*l2_ESO_d*b0*delta_t^2;
            -1/2*l3_ESO_d*b0*delta_t^2];

DLADRC_max_effort = 100-motor_bias;


%%
% %%%%%%%%%%%%%%% D-LQR-MADRC %%%%%%%%%%%%%%% %


% Getting a plant model in controllable canonical form
% From transfer function at first in continous time

Q = diag([1/(1)^2 1/(1.2)^2]);
R = 1/(8)^2;

k_ESO = 5;


A_m_lqr = [0 1;
           0 -a1];
b_m_lqr = [0;
           b0];

c_T_m_lqr = [1 0];

% Disturbance model
A_d_m_lqr = 0;
c_T_d_m_lqr = 1;

e_m = [0;
       1];

% Augmented Linear ESO (MESO)

A_o_lqr = [A_m_lqr e_m*c_T_d_m_lqr;
            0 0 A_d_m_lqr];
b_o_lqr = [b_m_lqr;
             0];
c_o_T_lqr = [c_T_m_lqr 0];

% Calculate Controller Gains
k_c_T_lqr = lqr(A_m_lqr,b_m_lqr,Q,R);
k_d_T_lqr = (e_m*c_T_d_m_lqr)\b_m_lqr;


% For Steady-State Performance
k_s_lqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_lqr)*b_m_lqr);


% C. State-Space System
countinous_ss_lqr = ss(A_m_lqr,b_m_lqr,c_T_m_lqr,0);

% Discretize State-Space and update controller gains
discrete_ss_lqr = c2d(countinous_ss_lqr,delta_t,'zoh');

A_m_dlqr = discrete_ss_lqr.A;
b_m_dlqr = discrete_ss_lqr.B;
c_m_dlqr = discrete_ss_lqr.C;
k_c_T_dlqr = dlqr(A_m_dlqr,b_m_dlqr,Q,R);
k_s_dlqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_dlqr)*b_m_lqr);

k_d_T_dlqr = k_d_T_lqr;
%k_d_T_dlqr = k_c_T_lqr;

% Discretize Observer
countinous_o_lqr = ss(A_o_lqr,b_o_lqr,c_o_T_lqr,0);
discrete_o_lqr = c2d(countinous_o_lqr,delta_t,'zoh');

A_o_d = discrete_o_lqr.A;
B_o_d = discrete_o_lqr.B;
C_o_d = discrete_o_lqr.C;


% Observer Gains
T_settle_98 = 1/k_ESO; % Seconds
N = 2;
omega_o = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

% z-domain pole locations
z_MESO = exp(-omega_o * delta_t);
l_MESO_d = acker(A_o_d',C_o_d',[z_MESO,z_MESO,z_MESO])';

A_MESO_d = A_o_d - l_MESO_d*C_o_d*A_o_d;
B_MESO_d = B_o_d - l_MESO_d*C_o_d*B_o_d;

DLQRMADRC_max_effort = 100-motor_bias;

%% Plain Discrete LQR For Comparison

% Controller Gains
k_dlqr = k_c_T_dlqr;

%k_s_lqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_lqr)*b_m_lqr);

% Observer Gains
T_settle_98 = 1/k_ESO; % Seconds
N = 2;
omega_o = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

% z-domain pole locations
z_o = exp(-omega_o * delta_t);
l_o_d = acker(A_m_dlqr',c_m_dlqr',[z_o,z_o])';


A_O_d = A_m_dlqr - l_o_d*c_m_dlqr*A_m_dlqr;
B_O_d = b_m_dlqr - l_o_d*c_m_dlqr*b_m_dlqr;

DLQR_max_effort = 100-motor_bias;

%% Discrete LQR-ADRC (No model) 

% Model in Continous Time (Integrator Chain)
A_lqr_adrcnm = [0 1;
         0 0];
B_lqr_adrcnm = [0;
         1];
C_lqr_adrcnm = [1 0];

% Discretize
integrator_chain_lqr_adrcnm = ss(A_lqr_adrcnm, B_lqr_adrcnm, C_lqr_adrcnm,0);
dis_integrator_chain_lqr_adrcnm = c2d(integrator_chain_lqr_adrcnm,delta_t,'zoh');

A_dlqr_adrcnm = dis_integrator_chain_lqr_adrcnm.A;
B_dlqr_adrcnm = dis_integrator_chain_lqr_adrcnm.B;

k_dlqr_adrcnm = dlqr(A_dlqr_adrcnm,B_dlqr_adrcnm,Q,R);
k_s_dlqr_adrcnm = -inv(C_lqr_adrcnm*inv(A_lqr_adrcnm-B_lqr_adrcnm*k_dlqr_adrcnm)*B_lqr_adrcnm);

k_d_T_dlqr_adrcnm = k_d_T_lqr; 

DLQR_ADRCNM_max_effort = 100-motor_bias;


%% Old Tunings

% Tunings For runs 20

%% %% Experimental %% %%

%% Lightweight Kalman

A_kalman = [1 delta_t;
            0    1  ]; 
C_kalman = eye(2);

Q_kalman = diag([0.000007, (0.07)^2]);
R_kalman = diag([(calc_sigma(0.8))^2, (0.07)^2]);

%% For Eval
inbuilt_reset_time = amfitrack_time + 1;
