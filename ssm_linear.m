function [X,fvec,X_spec,THD] = ssm_linear(a,x0,E,Mmt,Cmt,Rmt,...
                                            fs,fsig,inp_p,Rl,s_ind,e_ind)

% This function calculates the output voltage for given input signal and load
% resistance.
%
% Input Parameters
% a     -  Radius of the membrane
% x0   -  Equilibrium spacing between membrane and back plate
% E     -  Polarization voltage
% Mmt -  Total Mass
% Cmt  -  Total Compliance
% Rmt  -  Total Damping
% fs     -  Sampling frequency (Hz)
% fsig   -  Input signal frequency
% inp_p -  Discrete input signal
% Rl     -  Load resistance
% s_ind -  Starting index for the FFT computation
% e_ind -  Ending index for the FFT computation
%
% Output Parameters
% X      -   Matrix containing the time domain responses in each row in the
%             below order
%            1 - Charge
%            2 - Displacement
%            3 - Velocity
%            4 - Time varying Capacitance
%            5 - Output voltage
%            6 - Current
% fvec   -  Vector containing the frequencies of FFT
% X_spec - Matrix containing the FFT of variables under X
% THD    - Vector containing the THDs of variables under X

% constants
e0 = 8.85e-12;      % Faraday's Constant
S  = pi*a^2;        % Area of diaphragm
Ce0 = e0*S/x0;      % capacitence between diaphragm and backplate

% Initial values of state variables
q = 0;
x = 0;
u = 0;
X = [q;x;u]; % State vector

% F and G matrices
F = [-1/(Rl*Ce0)     -E/(Rl*x0)        0;
         0                0            1;
    -E/(x0*Mmt)       -1/(Mmt*Cmt)   -Rmt/Mmt];
G = [0;0;-S/Mmt];

% Compute the state variables
N = length(inp_p);
for ii = 1:N-1
    X(:,ii+1) = X(:,ii) + (1/fs)*(F*X(:,ii)+G*inp_p(ii));
end

% Find the time varying capacitence
X(4,:) = e0*S./(x0+X(2,:));

% Find the Output voltage
X(5,:) = E*X(2,:)/x0 + X(1,:)/Ce0;

% Find the current from Output voltage
X(6,:) = -X(5,:)/Rl;

% Charge
spec_t = fft(X(1,s_ind+1:e_ind));
Nfvec = length(spec_t);
fvec = (0:Nfvec/2-1)*fs/Nfvec;
X_spec(1,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% Displacement
spec_t = fft(X(2,s_ind+1:e_ind));
X_spec(2,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% Velocity
spec_t = fft(X(3,s_ind+1:e_ind));
X_spec(3,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% Time varying Capacitance
spec_t = fft(X(4,s_ind+1:e_ind));
X_spec(4,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% Output voltage
spec_t = fft(X(5,s_ind+1:e_ind));
X_spec(5,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% Current
spec_t = fft(X(6,s_ind+1:e_ind));
X_spec(6,:) = 2*(spec_t(1:Nfvec/2)/Nfvec);

% THD
nmax = 10; % No of harmonics to consider in the THD calculation
Nsl = 1;   % No of side lobes around the main lobe to consider 
THD = zeros(1,size(X,1));
for ii = 1:size(X,1)
    [THD(ii),~] = thd_nu(abs(X_spec(ii,:)),fsig,fvec,nmax,Nsl);
end

end