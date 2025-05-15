%LADRC: Linear ADRC in continuous time
clear
s = tf("s");

ex1_sys = 1 / (s * (s - 3));
ex2_sys = 3 / (s^2 + 6*s + 3);

K = 3;
D = 3;
T = 2;
test_plant = K / (T^2*s^2 + 2*T*D*s + 1); 

plant = test_plant;


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

T_settle_98 = 0.1; % Second

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


%% Filtering
ss_a0 = 1;
ss_a1 = 2*T*D;
ss_a2 = T^2;
ss_b0 = K;

ss_A = [     0      1      0;
             0      0      1;
        -ss_a0 -ss_a1 -ss_a2];

ss_B = [0;
        0;
      ss_b0];

ss_C = [1 0 0];

ss_D = 0;
sim_noise_pwr = 0.1;
kf_R = [sim_noise_pwr];
kf_Q = 10^-3 * eye(3);

test_ss = ss(test_plant)