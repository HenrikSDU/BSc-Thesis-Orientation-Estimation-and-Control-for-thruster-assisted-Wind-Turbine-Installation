%%
clc
clear

%% Functions
function R = Rzxy(phi,theta,psi)
% R = Rzxy(phi,theta,psi) computes the Euler angle
% rotation matrix R in SO(3) using the zyx convention
%
% Author:   Zhengru Ren
% Date:     18.9.2017



R = [   cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta),     -cos(phi)*sin(psi),     cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi);
        cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta),     cos(phi)*cos(psi),    	sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi);
    	-cos(phi)*sin(theta),                                   sin(phi),           	cos(phi)*cos(theta)];
end

%%
totalMass_s = 0.391;

%I_s = eye(3) .*[22.54e-3; 0.01995e-3; 22.54e-3];%::::::::::::::::::::::::::
I_s = eye(3) .*[0.0159; 0.00041299325; 0.0146]

M_s=[totalMass_s*eye(3) zeros(3); zeros(3) I_s];
inv_M_s =inv(M_s);

posCOG_s = [0;0.22;0];
posCOG_s = [0;0;0];

tot_liftWire = 0.7;%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::


L = 0.8;
W=0.05; %Width for approximative surface area

meanWSP = 1;

%% Coordinate
g = 9.81;


%% Crane tip
posCraneTip_s = [0;0;-1.2];%:::::::::::::::::::::::::::::::::::::::::::::::::
installation_height = posCraneTip_s+[0;0;tot_liftWire];

%% Lift wire and sling stiffness
E = 35e9;
D = 0.002;
A = pi/4*D^2;
EA = E*A;

gamma = 0.1;

lw1_s_init = tot_liftWire-0.14;%::::::::::::::::::::::::::::::::::::::::::::::::::::::::
k_w1_s = gamma * EA/lw1_s_init;

liftSpeed = 0;
%% Hook
m_h_s = 0.01;

elong_mainlw_s = (totalMass_s+m_h_s)*g/k_w1_s;

p_h_s_init = posCraneTip_s + [0;0;lw1_s_init + elong_mainlw_s];

%only for test in different init cond (disregardedelongation)
%p_h_s_init = R_b2n*(posCraneTip_s + [0;0;lw1_s_init]) 

p_h_dot_s_init = [0;0;0];

%% Blade
pitch=0;
yaw = 0;

p_s_init_n = [0;0;-0.51];%:::::::::::::::::::::::::::::::::::::::::
Theta_s_init_n = [0;pitch*pi/180;yaw*pi/180];

eta_s_init_n = [p_s_init_n;Theta_s_init_n];
nu_s_init_b = zeros(6,1);

R_b2n = Rzxy(Theta_s_init_n(1),Theta_s_init_n(2),Theta_s_init_n(3));

%% Slings

R_slings = [cosd(pitch) -sind(pitch); sind(pitch) cosd(pitch)];

slings_p1=[0.0425; -0.0525];
slings_p2=[-0.0425; -0.0525];

slings_p1_new = R_slings*slings_p1;
slings_p2_new = R_slings*slings_p2;%??????????????????????????why-to have same sling position when blade points down...



posSling_s_b = [posCOG_s + [slings_p1_new(1); -0.05; slings_p1_new(2)], ...
              posCOG_s + [slings_p1_new(1); 0.05; slings_p1_new(2)],...
              posCOG_s + [slings_p2_new(1); -0.05; slings_p2_new(2)],...
              posCOG_s + [slings_p2_new(1); 0.05; slings_p2_new(2)]];

for i =1:length(posSling_s_b(1,:))
    posSling_s_n(:,i) = eta_s_init_n(1:3)+R_b2n*(posSling_s_b(:,i)-posCOG_s);
end
posSling_s_n(:,:);  

k_w2_s = k_w1_s;
k_w3_s = k_w2_s;

elong_lw2_s = 0.5*(totalMass_s)*g/(sqrt(3)*k_w2_s); %%Modified to include 0.5 to account for 4 slings
elong_lw3_s = elong_lw2_s;

lw2_s = norm(p_h_s_init- posSling_s_n(:,1))-elong_lw2_s
lw3_s = norm(p_h_s_init- posSling_s_n(:,2))-elong_lw3_s;


d_w1_s = k_w1_s/10e3;
d_w2_s = k_w2_s/10e3;
d_w3_s = d_w2_s;

%% Thrusters

posThruster_s_b = [posCOG_s + [0.03; 0.21; -0.05], ... %+Yaw
              posCOG_s + [0.03; -0.21; -0.05],...      %-Yaw
              posCOG_s + [0; 0.1825; -0.05],...        %+Roll
              posCOG_s + [0; -0.1825; -0.05]];         %-Roll

for i =1:length(posThruster_s_b(1,:))
    posThruster_s_n(:,i) = eta_s_init_n(1:3)+R_b2n*(posThruster_s_b(:,i)-posCOG_s);
end
posThruster_s_n(:,:);
%% Simulation
tend = 1000;
stepTime = 100;

%% Control Vars
%%
% %%%%%%%%%%%%%%% LADRC %%%%%%%%%%%%%%% %
% Controller parameters

N = 2;
b0 =  3; 

% State space matrices of virtual plant
A_adrcl = [zeros(N, 1)  eye(N, N) ;
            0       zeros(1, N)];
b_adrcl = [zeros(N-1, 1);
             b0      ;
              1      ];

c_t = [1, zeros(1, N)];
%c = [1, 0, 0];
% Controller tuning based on desired time domain characteristics

T_settle_98 = 3; % Second

omega_cl = (1.9 + 2.09*N - 0.07*N^2) / T_settle_98;

k = ones(N,1); % controller_gains

for i=1:1:N

    k(i) = factorial(N) / (factorial(N-i + 1) * factorial(i-1)) * omega_cl^(N-i+1);

end 

% Tuning of Observer

K_ESO = 10;

omega_O = K_ESO * omega_cl;

l = ones(N+1, 1); % Observer gains

for i = 1:1:N+1
    
    l(i) = factorial(N+1) / (factorial(N-i+1)*factorial(i)) * omega_O^i;

end 

LADRC_max_effort = 30;
%%
% %%%%%%%%%%%%%%% DLADRC %%%%%%%%%%%%%%% %

b0 =  1.025; % Critical Gain Obtained in ParameterIdentification.m
delta_t = 0.01;

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

DLADRC_max_effort = 30;

%%
% %%%%%%%%%%%%%%% D-LQR-MADRC %%%%%%%%%%%%%%% %

delta_t = 0.01;
a1 = 0.16072; % Drag Coefficient Estimate obtained in ParameterIdentification.m
b0 =  1.025; % Critical Gain Obtained in ParameterIdentification.m
% Getting a plant model in controllable canonical form
% From transfer function at first in continous time

Q = diag([1/(3)^2 1/(0.9)^2]);
R = 1/(2)^2;

k_ESO = 7;


A_m_lqr = [0 1;
         -a1 0];
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
k_d_T_lqr = (e_m*c_T_d_m_lqr) \ b_m_lqr;


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
k_d_T_dlqr = k_c_T_lqr;

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

DLQRMADRC_max_effort = 30;