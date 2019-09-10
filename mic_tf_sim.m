function [e_ref,e_noref,ph_ref,ph_noref,PdB,Ts] = mic_tf_sim(Pi,f)

% This program uses the measured resonant frequency and sensitivity from
% LAB B of 31220 course and computes the theoretical lumped parameter
% values for B&K 4133 microphone.
%
% Input Parameters
% Pi - Level of Input Signal
% f  - Frequency vector
% 
% Output Parameters
% e_ref    - Magnitude response in complex form including the membrane reflections.
% e_noref  - Magnitude response in complex form excluding the membrane reflections.
% ph_ref   - Phase response in degrees including the membrane reflections.
% ph_noref - Phase response in degrees excluding the membrane reflections.
% Pdb      - Sound Pressure Level in dB re 20uPa
% Ts       - Transfer function simulating the membrane reflections in
%            complex form.

a = 8.95e-3/2;         % Radius of membrane
S = pi*a^2;            % Area of membrane
x0 = 20.77e-6;         % equlibrium distance between membrane and backplate
Cmd = 20*1e-6;         % Mechanical compliance of membrane
Mmd = 1.5*1e-6;        % Mechanical mass of the membrane
E = 200;               % Polarization voltage
rho = 1.2;             % Density of air
c = 344;               % Speed of sound

% Adjusted values from LTspice to match measured response
Cab2 = 2*1e-12;
Mas = 600;
Ras = 400*1e+6;

% We have mechanical and acoustical compliances so let's find total
% compliance of the system first
Cmt = 1/((1/Cmd)+(S^2/Cab2));

% Compute the radiation impedance components
Ma1 = 8*rho/(3*pi^2*a);
Ra1 = 0.441*rho*c/(pi*a^2);
Ra2 = (rho*c)/(pi*a^2);
Ca1 = 5.94*a^3/(rho*c^2);

% Compute the coeffecients of T(s)
Rpar1 = Ra1*Ra2/(Ra1+Ra2);
Rpar2 = (Ra1+Ra2)*Ra2/(Ra1+2*Ra2);
b1 = Rpar1*Ca1 + Ma1/Rpar2 ;
b2 = 2*Ra1*Ma1*Ca1/(Ra1+Ra2) ;
c1 = Rpar1*Ca1 + Ma1/(Ra1+Ra2) ;
c2 = b2/2 ;

% Compute T(s)
w = 2*pi*f;
s = 1j*w;
Ts = (1+b1*s+b2*s.^2)./(1+c1*s+c2*s.^2);

% Compute Mmt and Rmt
Mmt = Mmd + (Mas+Ma1)*S^2;
Rmt = S^2*Ras;

% Compute fres and Q
fres = 1/(2*pi*sqrt(Mmt*Cmt));
Q = (2*pi*fres*Mmt)/Rmt;

% Compute the Transfer function
ws = 2*pi*fres;
lp_tf = 1./((s/ws).^2+(1/Q)*(s/ws)+1); % Low Pass TF part of actual TF
M = E*S*Cmt/x0;
e_ref = -M*Pi*(lp_tf.*Ts);
e_noref = -M*Pi*(lp_tf);

% Compute the Phase
ph_ref = atan2d(imag(e_ref),real(e_ref));
ph_noref = atan2d(imag(e_noref),real(e_noref));

% Send out the SPL
PdB = 20*log10(Pi/20e-6);
end