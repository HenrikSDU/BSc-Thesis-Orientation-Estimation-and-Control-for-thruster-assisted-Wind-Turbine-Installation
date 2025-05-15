clc
clear

%% Rotation Matrix
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

%% Wind field
meanWindSpeed = 12;%%Modified
load(['D:/AAStudium/Study/Semester_6/BScThesis/Marin_blade_v1_0/data/MannTurbExample/WSP10_ClassC_Seed48_1000s.mat']);

%% Coordinate
load('D:/AAStudium/Study/Semester_6/BScThesis/Marin_blade_v1_0/data/bladePar.mat');
g = 9.81;
installation_height = -90;
m_yoke = 20e3;

%% Crane tip
posCraneTip = [0;0;-110];

%% Lift wire and sling stiffness
E = 2.1e11;
D = 0.25;
A = pi/4*D^2;
EA = E*A;

gamma = 0.1;

lw1_init = 15; %% Original 10
k_w1 = gamma * EA/lw1_init;

liftSpeed = 0;
%% Hook
m_h = 1e3;

elong_mainlw = (bladePar.totalMass+m_h+m_yoke)*g/k_w1; %Hooks Law badum tss

p_h_init = posCraneTip + [0;0;lw1_init + elong_mainlw];
p_h_dot_init = [0;0;0];

%% Blade
pitch=0;
yaw = 0;

damping_en =1;
pitch_damping=100e+3;

d_sling=5; %Distance in -z_n direction of sling attachement point

p_init_n = [0;0;-90];
Theta_init_n = [0;pitch*pi/180;yaw*pi/180];%%%%%% Modified!

eta_init_n = [p_init_n;Theta_init_n];
nu_init_b = zeros(6,1); % Velocities eta_d

R_b2n = Rzxy(Theta_init_n(1),Theta_init_n(2),Theta_init_n(3));

%% Slings
posSling_b = [bladePar.posCOG + [-2;-4.5;0] , bladePar.posCOG + [-2;4.5;0]] 
%test slings_delta_xyz=[d_sling*sind(pitch);4.5;d_sling*cosd(pitch)]

% posSling_b = [bladePar.posCOG + [d_sling*sind(pitch);-4.5;-d_sling*cosd(pitch)], ...
%               bladePar.posCOG + [d_sling*sind(pitch);4.5;-d_sling*cosd(pitch)],...
%               bladePar.posCOG + [d_sling*sind(pitch);4.5;-d_sling*cosd(pitch)],...
%               bladePar.posCOG + [d_sling*sind(pitch);4.5;-d_sling*cosd(pitch)]]%%% Rotated to match horizontal blade rotation% original z=-2

R_slings = [cosd(pitch) -sind(pitch); sind(pitch) cosd(pitch)];

slings_p1=[3; -2];
slings_p2=[-3; -2];

slings_p1_new = R_slings*slings_p1
slings_p2_new = R_slings*slings_p2



posSling_b = [bladePar.posCOG + [slings_p1_new(1); -4.5; slings_p1_new(2)], ...
               bladePar.posCOG + [slings_p1_new(1); 4.5; slings_p1_new(2)],...
               bladePar.posCOG + [slings_p2_new(1); -4.5; slings_p2_new(2)],...
               bladePar.posCOG + [slings_p2_new(1); 4.5; slings_p2_new(2)]]

% posSling_b = [bladePar.posCOG + [3; -4.5; -2], ...
%               bladePar.posCOG + [3; 4.5; -2],...
%               bladePar.posCOG + [-3; -4.5; -2],...
%               bladePar.posCOG + [-3; 4.5; -2]]%%% Rotated to match horizontal blade rotation% original z=-2



for i =1:length(posSling_b(1,:))
    posSling_n(:,i) = eta_init_n(1:3)+R_b2n*(posSling_b(:,i)-bladePar.posCOG);
end


k_w2 = 1e8;
k_w3 = k_w2;

elong_lw2 = 0.5*(bladePar.totalMass+m_yoke)*g/(sqrt(3)*k_w2); %%Modified to include 0.5 to account for 4 slings
elong_lw3 = elong_lw2;

lw2 = norm(p_h_init- posSling_n(:,1))-elong_lw2;
lw3 = norm(p_h_init- posSling_n(:,2))-elong_lw3;

%% Taglines
posTuggerline_b = [bladePar.posCOG+[0;-4.5;0],bladePar.posCOG+[0;-35;0]]; %Modified original: 0;+-4.5;0

for i =1:length(posTuggerline_b(1,:))
    posTuggerline_n(:,i) = eta_init_n(1:3)+R_b2n*(posTuggerline_b(:,i)-bladePar.posCOG);
    posTuggerlineBase_n(:,i) = posTuggerline_n(:,i) - [10;0;0];
end

lt1_init = norm(posTuggerline_n(:,1)- posTuggerlineBase_n(:,1));
lt2_init = norm(posTuggerline_n(:,2)- posTuggerlineBase_n(:,2));

k_t1 = 1e7;
k_t2 = k_t1;

for i =1:length(posSling_b(1,:))
    pLiftwire_init_n(:,i) = eta_init_n(1:3)+R_b2n*(posSling_b(:,i)-bladePar.posCOG);
end

d_w1 = k_w1/10e3;
d_w2 = k_w2/10e3;
d_w3 = d_w2;
d_t1 = 0; %k_t1/1e4;
d_t2 = d_t1;

%% Simulation
tend = 1000;
stepTime = 100;


%% %%%%% CONTROL %%%%% %%

%%
% %%%%%%%%%%%%%%% LADRC %%%%%%%%%%%%%%% %
% Controller parameters

N = 2;
b0 =  9.345442381739944e-04; 

% State space matrices of virtual plant
A_adrcl = [zeros(N, 1)  eye(N, N) ;
            0       zeros(1, N)];
b_adrcl = [zeros(N-1, 1);
             b0      ;
              1      ];

c_t = [1, zeros(1, N)];
%c = [1, 0, 0];
% Controller tuning based on desired time domain characteristics

T_settle_98 = 0.2; % Second

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

LADRC_max_effort = 900;

%% ADRC + Additional Model Information

% Getting a plant model in controllable canonical form
% From transfer function

a1 = 0.14453;

A_m = [0 1;
       0 -a1];
b_m = [0;
       b0];

c_T_m = [1 0];


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

%% LQR-ADRC + Additional Model Information

% Getting a plant model in controllable canonical form
% From transfer function
Q = diag([1/(0.001)^2 1/(0.1)^2]);
R = 1/(750)^2;

A_m_lqr = [0 1;
           0 -a1];
b_m_lqr = [0;
          b0];

c_T_m_lqr = [1 0];

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

%% LQR-ADRC with no model info
A_m_nmlqr = [0 1;
           0 0];
b_m_nmlqr = [0;
          b0];

c_T_m_nmlqr = [1 0];

% Disturbance model
A_d_m_nmlqr = 0;
c_T_d_m_nmlqr = 1;

e_m = [zeros(N-1,1);
            1     ];

% Augmented Linear ESO

A_o_nmlqr = [A_m_nmlqr e_m*c_T_d_m_nmlqr;
            0 0 A_d_m_nmlqr];
b_o_nmlqr = [b_m_nmlqr;
                0];
c_o_T_nmlqr = [c_T_m_nmlqr 0];


k_c_T_nmlqr = lqr(A_m_nmlqr,b_m_nmlqr,Q,R);
k_d_T_nmlqr = (e_m*c_T_d_m_nmlqr) \ b_m_nmlqr;

k_o_nmlqr = acker(A_o_nmlqr',c_o_T_nmlqr',[-omega_O,-omega_O,-omega_O])';

k_s_nmlqr = -inv(c_T_m_nmlqr*inv(A_m_nmlqr-b_m_nmlqr*k_c_T_nmlqr)*b_m_nmlqr);


%% Plain LQR based on simplified model

A_o_plqr = A_m_lqr;
b_o_plqr = b_m_lqr;
c_o_T_plqr = c_T_m_lqr;
% Controller Gains
k_plqr = lqr(A_m_lqr, b_m_lqr, Q, R);
k_s_plqr = -inv(c_T_m*inv(A_m-b_m*k_plqr)*b_m);
% Observer gains
k_o_plqr = acker(A_o_plqr',c_o_T_plqr',[-omega_O,-omega_O])';

%% Pole comparison
%countinous_controlled_ss_lqr = ss(A_m_lqr-k_c_T_lqr*b_m_lqr,b_m_lqr,c_T_m_lqr,0);
%pzmap(countinous_controlled_ss_lqr) % k_c_T_lqr
