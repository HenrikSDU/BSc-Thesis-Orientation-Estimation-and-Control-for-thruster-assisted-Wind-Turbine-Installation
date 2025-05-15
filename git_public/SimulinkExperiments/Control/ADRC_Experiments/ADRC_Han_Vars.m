clear; close all; clc;
%% Variables for ADRC alá Han


% Test Plant (2nd Order)
s = tf("s");
K = 10;
D = 5;
T = 3;
test_plant = K / (T^2*s^2 + 2*T*D*s + 1); 

plant = test_plant;


% General Params
h = 0.01; % Sample time
r = 100000; % Actuator Limit

% Controller Params
h0 = 8*h; % Controller Param = Sample time
r0 = 2; % Limits actuator / Acceleration Limit for trajectory generation
c = 50; % Fine Tuning Param for NLSEF
b0 =  K/T^2;  % Critical Gain: Formula only holds for second order system


% Observer Params
alpha_1 = 1;
alpha_2 = 1;

beta_01 = 1;
beta_02 = 1 / (2*h^0.5);
beta_03 = 2 / (5^2 * h^1.2);

delta = 0.05;

%% Linear ADRC for comparison


% Controller parameters

N = order(plant);
b0 =  K/T^2;  % Only holds for second order system

% State space matrices of virtual plant
A = [zeros(N, 1)  eye(N, N) ;
            0       zeros(1, N)];
b = [zeros(N-1, 1);
             b0      ;
              1      ];

c_t = [1, zeros(1, N)];
%c = [1, 0, 0];
% Controller tuning based on desired time domain characteristics

T_settle_98 = 1; % Second

omega_cl = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

k = ones(N,1); % controller_gains

for i=1:1:N

    k(i) = factorial(N) / (factorial(N-i + 1) * factorial(i-1)) * omega_cl^(N-i+1);

end 

% Tuning of Observer

K_ESO = 7;

omega_O = K_ESO * omega_cl;

l = ones(N+1, 1); % Observer gains

for i = 1:1:N+1
    
    l(i) = factorial(N+1) / (factorial(N-i+1)*factorial(i)) * omega_O^i;

end 

%% ADRC + Additional Model Information

% Getting a plant model in controllable canonical form
% From transfer function

A_m = [zeros(N-1,1) eye(N-1,N-1);
       -1/T^2 -2*T/(T^2)];
b_m = [zeros(N-1,1);
           b0     ];

c_T_m = [1 zeros(1,N-1)];

% Disturbance model
A_d_m = 0;
c_T_d_m = 1;

e_m = [zeros(N-1,1);
            1     ];

% Augmented Linear ESO

A_o = [A_m e_m*c_T_d_m;
            0 0 A_d_m];
b_o = [b_m;
        0];
c_o_T = [c_T_m 0];

k_c_T = acker(A_m,b_m, [-omega_cl,-omega_cl]);
k_d_T = (e_m*c_T_d_m) \ b_m;

k_o = acker(A_o',c_o_T',[-omega_O,-omega_O,-omega_O])';

k_s = -inv(c_T_m*inv(A_m-b_m*k_c_T)*b_m);


%% Discrete LADRC

delta_t = 0.01;

N = 2;
k_ESO = 7; 
T_settle_98 = 1; % Seconds
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

%% Plain LQR-ADRC 
A_lqr = [0 1;
         0 0];
B_lqr = [0;
         1];
C_lqr = [1 0];
Q = diag([1/(0.1)^2 1]);
R = 1;

K_lqr = lqr(A_lqr,B_lqr,Q,R);
K_s_lqr = -inv(C_lqr*inv(A_lqr-B_lqr*K_lqr)*B_lqr);



%% LQR-ADRC + Additional Model Information

% Getting a plant model in controllable canonical form
% From transfer function

A_m_lqr = [zeros(N-1,1) eye(N-1,N-1);
       -1/T^2 -2*T/(T^2)];
b_m_lqr = [zeros(N-1,1);
           b0     ];

c_T_m_lqr = [1 zeros(1,N-1)];

% Disturbance model
A_d_m_lqr = 0;
c_T_d_m_lqr = 1;

e_m = [zeros(N-1,1);
            1     ];

% Augmented Linear ESO

A_o_lqr = [A_m_lqr e_m*c_T_d_m_lqr;
            0 0 A_d_m_lqr];
b_o_lqr = [b_m_lqr;
        0];
c_o_T_lqr = [c_T_m_lqr 0];


k_c_T_lqr = lqr(A_m_lqr,b_m_lqr,Q,R);
k_d_T_lqr = (e_m*c_T_d_m_lqr) \ b_m_lqr;

k_o_lqr = acker(A_o_lqr',c_o_T_lqr',[-omega_O,-omega_O,-omega_O])';

k_s_lqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_lqr)*b_m_lqr);


% C. State Space System
countinous_ss_lqr = ss(A_m_lqr,b_m_lqr,c_T_m_lqr,0);

%figure;
%step(countinous_ss_lqr,10)

countinous_controlled_ss_lqr = ss(A_m_lqr-k_c_T_lqr*b_m_lqr,b_m_lqr,c_T_m_lqr,0);

%figure;
%step(countinous_controlled_ss_lqr,10)

% C. Observer

countinous_o_lqr = ss(A_o_lqr,b_o_lqr,c_o_T_lqr,0);

%% LQR-ADRC + Additional Model Information --> Discrete

discrete_ss_lqr = c2d(countinous_ss_lqr,delta_t,'zoh');

A_m_dlqr = discrete_ss_lqr.A;
b_m_dlqr = discrete_ss_lqr.B;
c_m_dlqr = discrete_ss_lqr.C;
k_c_T_dlqr = dlqr(A_m_dlqr,b_m_dlqr,Q,R);
k_d_T_dlqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_dlqr)*b_m_lqr);
%k_s_lqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_lqr)*b_m_lqr);

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


%% Plain Discrete LQR

% Discretize plant model
plain_ss = ss(A_m, b_m,c_T_m,0);
discrete_plain_ss = c2d(plain_ss,delta_t,'zoh');

A_m_d = discrete_plain_ss.A;
B_m_d = discrete_plain_ss.B;
C_m_d = discrete_plain_ss.C;

% Controller Gains
k_dlqr = dlqr(A_m_d, B_m_d, Q, R);
k_s_dlqr = -inv(c_T_m*inv(A_m-b_m*k_dlqr)*b_m);
%k_s_lqr = -inv(c_T_m_lqr*inv(A_m_lqr-b_m_lqr*k_c_T_lqr)*b_m_lqr);

% Observer Gains
T_settle_98 = 1/k_ESO; % Seconds
N = 2;
omega_o = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

% z-domain pole locations
z_o = exp(-omega_o * delta_t);
l_o_d = acker(A_m_d',C_m_d',[z_o,z_o])';


A_O_d = A_m_d - l_o_d*C_m_d*A_m_d;
B_O_d = B_m_d - l_o_d*C_m_d*B_m_d;


%% Discrete simul
test_ss = ss(A_m,b_m,c_T_m,0);
d_test_ss = c2d(test_ss,0.01);
A_d_test = d_test_ss.A; 
B_d_test = d_test_ss.B;
C_d_test = d_test_ss.C;
