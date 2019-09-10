function eig_val = ssm_linear_stab_check(a,x0,E,Mmt,Cmt,Rmt,stabf,Rl)

% This function calculates the eigen values for given frequencies and load
% resistance which are used to evaluate the stability of a state space
% model
%
% Input Parameters
% a     -  Radius of the membrane
% x0   -  Equilibrium spacing between membrane and back plate
% E     -  Polarization voltage
% Mmt -  Total Mass
% Cmt  -  Total Compliance
% Rmt  -  Total Damping
% Stabf -  Vector of frequencies to evaluate the stability at
% Rl     -  Vector of load resistances to evaluate the stability at
%
% Output Parameters
% eig_val - Matrix with each row having the absolute of the maximum eigen
%             value for all frequencies of a load resistance.

% Constants
e0 = 8.85e-12;      % Faraday's Constant
S  = pi*a^2;         % Area of diaphragm
Ce0 = e0*S/x0;      % capacitence between diaphragm and damping material

% Initialize the variables depending on the size of inputs
Nf = length(stabf);
NR = length(Rl);
eig_val = zeros(NR,Nf);
I = eye(3);

% Loop through each value of load resistance
for ii = 1:NR
    
    % F Matrix
    F = [-1/(Rl(ii)*Ce0)     -E/(Rl(ii)*x0)        0;
               0                    0              1;
          -E/(x0*Mmt)       -1/(Mmt*Cmt)   -Rmt/Mmt];
    
    % Loop through each frequency
    for jj = 1:Nf
        eig_mat = I+F/stabf(jj) ;
        eig_val(ii,jj) = max(abs(eig(eig_mat)));
    end

end

end