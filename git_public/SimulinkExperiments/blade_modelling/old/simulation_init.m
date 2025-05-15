%% Initialization for simulation
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

%% Wind field and Wind Loads
A_b_s = 0.040; % Approximate area of blade
C_d_s = 1.050; % Drag coefficient
rho_s = 1.1225; % Air Density
meanWSP = 2;
% "D:\AAStudium\Study\Semester_6\BScThesis\Marin_blade_v1_0\data\MannTurbExample\WSP10_ClassC_Seed48_1000s.mat"
% load(['D:/AAStudium/Study/Semester_6/BScThesis/Marin_blade_v1_0/data/MannTurbExample/WSP10_ClassC_Seed48_1000s.mat']);

%% Model Parameters
totalMass_s = 0.391;
I_s = eye(3) .* [0.0159; 0.0004199235; 0.0146];

M_s = [totalMass_s*eye(3) zeros(3); zeros(3) I_s];
inv_M_s = inv(M_s);

posCOG_s = [0; 0; 0];
tot_liftWire = 0.7;
L = 0.8;
W = 0.05; % Width for approximate surface area



%% Constants, Coordinates and Tower Crane Parameters
% load('D:/AAStudium/Study/Semester_6/BScThesis/Marin_blade_v1_0/data/bladePar.mat');
g = 9.81;
posCraneTip_s = [0; 0; -1.2];
installation_height = posCraneTip_s + [0; 0; tot_liftWire];


%% Lift wire and sling stiffness
E = 35e9;
D = 0.002;
A = pi/4*D^2;
EA = E*A;

gamma = 0.1;

lw1_s_init = tot_liftWire - 0.14;
k_w1_s = gamma * EA/lw1_s_init;

liftSpeed = 0;
%% Hook
m_h_s = 0.01;

elong_mainlw_s = (totalMass_s + m_h_s) * g / k_w1_s;

p_h_s_init = posCraneTip_s + [0;0;lw1_s_init + elong_mainlw_s];
p_h_dot_s_init = [0;0;0];

%% Blade
pitch = 0;
yaw = 0;

p_s_init_n = [0;0;-0.51];
Theta_s_init_n = [0; pitch * pi/180; yaw * pi/180];

eta_s_init_n = [p_s_init_n;Theta_s_init_n];
nu_init_b = zeros(6,1);

R_b2n = Rzxy(Theta_s_init_n(1),Theta_s_init_n(2),Theta_s_init_n(3));

%% Slings

R_slings = [cosd(pitch) -sind(pitch); 
            sind(pitch)  cosd(pitch)];

slings_p1 = [0.0425; -0.0525];
slings_p2 = [-0.0425; -0.0525];

slings_p1_new = R_slings * slings_p1;
slings_p2_new = R_slings * slings_p2;

posSling_s_b = [posCOG_s + [slings_p1_new(1); -0.05; slings_p1_new(2)], ...
                posCOG_s + [slings_p1_new(1); 0.05; slings_p1_new(2)], ...
                posCOG_s + [slings_p2_new(1); -0.05; slings_p2_new(2)], ...
                posCOG_s + [slings_p2_new(1); 0.05; slings_p2_new(2)]];


for i =1:length(posSling_s_b(1,:))
    posSling_s_n(:,i) = eta_s_init_n(1:3)+R_b2n*(posSling_s_b(:,i) - posCOG_s);
end



k_w2_s = k_w1_s;
k_w3_s = k_w2_s;

k_w_s = k_w1_s;
d_w_s = 10 * 10^3; % Dampening of Wire 

elong_lw2_s = (0.5*(totalMass_s))*g/(sqrt(3)*k_w2_s);
elong_lw3_s = elong_lw2_s;

lw2_s = norm(p_h_s_init- posSling_s_n(:,1)) - elong_lw2_s;
lw3_s = norm(p_h_s_init- posSling_s_n(:,2)) - elong_lw3_s;

%% Thrusters

posThrusters_s_b = [posCOG_s + [0.03;    0.21; -0.05], ... % -Yaw
                    posCOG_s + [0.03;   -0.21; -0.05], ... % +Yaw
                    posCOG_s + [   0;  0.1825; -0.05], ... % -Roll
                    posCOG_s + [   0; -0.1825; -0.05],     % +Roll
                    ];

for i = 1:length(posThrusters_s_b(1,:))
    posThrusters_s_n(:,i) = eta_s_init_n(1:3) +R_b2n * (posThrusters_s_b (:,i) - posCOG_s);
end
posThrusters_s_n(:,:);


max_thrust_g = 8;
max_thrust_N = 8 / 1000 * 9.81;

%% Simulation
tend = 1000;
stepTime = 100;