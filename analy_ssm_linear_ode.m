clc
clear
close all

% All the microphone parameters are actually controlled inside the function
% linear_ode45.m

%% Get the output from ODE45

fs = 96000;    % Sampling frequency
A  = 1;         % Signal level
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
inp_p = A*sin(2*pi*fsig*t);
s_ind = 1*fs;
e_ind = T*fs;
Rl = 10e+9;
N = T*fs;
x0 = 20.77e-6;
a = 8.95e-3/2;
E = 200;

tic
[X1,fvec1,X1_spec,THD1] = ssm_linear_ode(a,x0,E,t,fs,fsig,N,Rl,s_ind,e_ind);
eta1 = toc;

%% Get the output using Forward Euler

a = 8.95e-3/2;           % Radius of diaphragm
x0 = 20.77e-6;          % Equilibrium separation b/w plates when E = 0
E = 200;                  % Polarization voltage
Mmt = 2.09e-6;         % Total Mass
Cmt = 1.92e-5;         % Total Compliance
Rmt = 1.08;             % Total Damping
e0 = 8.85e-12;      % Faraday's Constant
S  = pi*a^2;        % Area of diaphragm
Ce0 = e0*S/x0;      % capacitence between diaphragm and backplate

fs = 231500;    % Sampling frequency
A  = 1;           % Signal level
AdB = 20*log10(A/20e-6);  % Signal level in dB
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*fsig*t);
s_ind = 1*fs;
e_ind = T*fs;
Rl = 10e+9;  % Load resistance

tic
[X2,fvec2,X2_spec,THD2] = ssm_linear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,p_in,Rl,s_ind,e_ind);
eta2 = toc;
%% Plot and Compare the Spectra from both methods

figure
semilogx(fvec1,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec2,20*log10(abs(X2_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+"dB re 20\muPa)")
legend('ODE45','Forward Euler')
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])