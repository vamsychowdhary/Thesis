function dXdt = linear_ode45(t,y)

% This function is called from ssm_linear_ode
%
% Input Parameters
% t - 2 row Vector containing current and next samples of time
% y - 3 column vector containing the current values of states
%
% Output Parameters
% dXdt - Differential change between states

% Microphone paraeters
a = 8.95e-3/2;           % Radius of diaphragm
x0 = 20.77e-6;          % Equilibrium separation b/w plates when E = 0
E = 200;                  % Polarization voltage
Mmt = 2.09e-6;         % Total Mass
Cmt = 1.92e-5;         % Total Compliance
Rmt = 1.08;             % Total Damping
e0 = 8.85e-12;      % Faraday's Constant
S  = pi*a^2;        % Area of diaphragm
Ce0 = e0*S/x0;      % capacitence between diaphragm and backplate
A = 1;
fsig = 20;
inp_p = A*sin(2*pi*fsig*t);
Rl = 10e+9;

% Compute the vector
dXdt = [-(y(1)/(Rl*Ce0))-(E*y(2)/(Rl*x0));
        y(3);
       -(E*y(1)/(x0*Mmt))-(y(2)/(Mmt*Cmt))-(y(3)*Rmt/Mmt)-(S*inp_p/Mmt)];
end