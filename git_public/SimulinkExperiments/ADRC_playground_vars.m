%%%% ADRC Playground
%% Constants

% Plant Constants
a1 = 2;
a2 = 3;
b = 0.2;

% Controller Parameters
b0 = 0.4;
h1 = 0.01;
c = 1;

r = 5;

% Sampling time 
h = 0.01;

beta1 = 1;
beta2 = 1 / (2*h^0.5);
beta3 = 2 / (5^2 * h^1.2);