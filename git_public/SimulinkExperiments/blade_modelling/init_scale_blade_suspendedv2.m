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

meanWSP = 1.6;

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
k_w1_s = 10 * gamma * EA/lw1_s_init;

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


d_w1_s =  10e3;%k_w1_s/10e3;
d_w2_s =  10e3;%k_w2_s/10e3;
d_w3_s = d_w2_s;

%% thrusters

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

