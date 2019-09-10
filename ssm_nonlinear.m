function [X,fvec,X_spec,THD_TIMD] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                                   fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL)

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
% f1     -  Modulation frequency for TIMD
% f2     -  Fundamental frequnecy for TIMD
% p_in  -  Discrete input signal
% mode - "THD" or "TIMD" depending on the input signal
% Rl     -  Load resistance
% s_ind -  Starting index for the FFT computation
% e_ind -  Ending index for the FFT computation
% qNL   -  Consider the nonlinearity under q
% uNL   -  Consider the nonlinearity under u
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
% THD_TIMD    - Vector containing the THDs of variables under X

% constants
e0 = 8.85e-12;      % Faraday's Constant
S  = pi*a^2;        % Area of diaphragm
Ce0 = e0*S/x0;      % capacitence between diaphragm and backplate

% Initial values of state variables
q = 0;
x = 0;
u = 0;
X = [q;x;u]; % State vector

% Convolve with the impulse response of the diaphragm reflections
N = length(p_in);

% G Matrix
G = [0;0;-S/Mmt];

% Compute the state variables
for ii = 1:N-1
    
    % Check for qNL
    if (qNL == "Y")
        F11 = -(1/(Rl*Ce0))*(1 + X(2,ii)/x0);
    else
        F11 = -1/(Rl*Ce0);
    end
    
    
    % Check for uNL
    if (uNL == "Y")
        F31 = -(1/(x0*Mmt))*(E + X(1,ii)/(2*Ce0));
    else
        F31 = -E/(x0*Mmt);
    end    
    
    % Update F matrix each time
    F = [F11     -E/(Rl*x0)                0;
            0              0                     1;
           F31     -1/(Mmt*Cmt)    -Rmt/Mmt];   
    
    X(:,ii+1) = X(:,ii) + (1/fs)*(F*X(:,ii)+G*(p_in(ii)));
    
end

% Find the time varying capacitence
X(4,:) = e0*S./(x0+X(2,:));

% Find the Output voltage
if (qNL == "Y")
    X(5,:) = E*X(2,:)/x0 + (X(1,:)/Ce0).*(1+X(2,:)/x0);
else
    X(5,:) = E*X(2,:)/x0 + X(1,:)/Ce0;  
end

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

if (mode == "THD")
    for ii = 1:6
        [THD_TIMD(ii),~] = thd_nu(abs(X_spec(ii,:)),fsig,fvec,nmax,Nsl);
    end
else
    for ii = 1:6
            [THD_TIMD(ii),~] = timd_nu(abs(X_spec(ii,:)),f1,f2,fs,fvec,nmax,Nsl);
    end   
end

end