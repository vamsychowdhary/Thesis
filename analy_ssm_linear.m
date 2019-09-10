% In this script, the transient analysis using a linear state space model
% is carried out to validate the model and condenser microphone behaviour

clc
clear
close all

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

%% Transient analysis

fs = 800000;    % Sampling frequency
A  = 1;           % Signal level
AdB = 20*log10(A/20e-6);  % Signal level in dB
fsig = 8000;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*fsig*t);
s_ind = 1*fs;
e_ind = T*fs;
Rl = 10e+9;  % Load resistance

[X,fvec,X_spec,THD] = ssm_linear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,p_in,Rl,s_ind,e_ind);

%% Plot the state variables
s_ind = 1*fs;
e_ind = 1.5*fs;
figure
subplot(2,1,1)
yyaxis left
plot(t(s_ind+1:e_ind),X(2,s_ind+1:e_ind)*1e+6,'LineWidth',1.2)
ylabel('Displacement (\mum)')
hold on
yyaxis right
plot(t(s_ind+1:e_ind),X(3,s_ind+1:e_ind)*1e+6,'LineWidth',1.2)
ylabel('Velocity (\mum/s)')
xlabel('Time (seconds)')
title("State Variables (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
set(gca,'FontSize',12)

subplot(2,1,2)
plot(t(s_ind+1:e_ind),X(1,s_ind+1:e_ind),'LineWidth',1.2,'color','r')
ylabel('Charge (C)')
xlabel('Time (seconds)')
title("State Variables (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Get the output voltage from MATLAB TF Simulation

% Get the frequency response with the same number of points of that LTSpice
% output
P_in = 1;
[~,e_noref,~,~,~,~] = mic_tf_sim(P_in,fvec);

% Apply ifft to get Impulse Response
E_noref = [e_noref fliplr(e_noref)];
mic_ir = ifft(E_noref,'symmetric'); 
fsn = 2*ceil(max(fvec));
tr = 0:1/fsn:2-1/fsn;
f_ir = 0:fsn/2 - 1;

% Convolution of IR with input signal
e_tf_amp = conv(p_in,mic_ir);     % Transient response of the TF output
e_tf_freq = fft(e_tf_amp(fsn+1:2*fsn));  % Only take till the length of input signal
Y = abs(e_tf_freq(1:fsn/2));
nmax = 10;
Nsl = 1;
[THD_tf,~] = thd_nu(Y,fsig,f_ir,nmax,Nsl);

% Get the time scale for the MATLAB simulated output
T = length(e_tf_amp)/fsn;
tn = 0:1/fsn:T-1/fsn;

%% Compare State Space Model and MATLAB IR Simulation
s_ind = 1*fs;
e_ind = 1.5*fs;
maxssm = max(X(5,s_ind+1:e_ind)*1e+3);
figure
subplot(3,1,1)
plot(t(s_ind+1:e_ind),X(5,s_ind+1:e_ind)*1e+3,'Color','b','LineWidth',1)
hold on
s_ind = find(tn >= 1,1);
e_ind = find(tn >= 1.5,1);
maxtf = max(e_tf_amp(s_ind:e_ind)*1e+3);
plot(tn(s_ind:e_ind),e_tf_amp(s_ind:e_ind)*1e+3,'--','Color','r','LineWidth',1)
ylabel('Voltage (mV)')
xlabel('Time (seconds)')
title("Output Voltage Comparison (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
legend('State Space Model','MATLAB IR Simulation')
set(gca,'FontSize',12)
%set(gcf,'Position',[100 100 900 250])
diff = maxtf - maxssm;

s_ind = 1*fs;
e_ind = 1.5*fs;
subplot(3,1,2)
yyaxis left
plot(t(s_ind+1:e_ind),(E*X(2,s_ind+1:e_ind)/x0)*1e+3,'LineWidth',1.2)
ylabel('Ex/x_0 term (mV)')
hold on
yyaxis right
plot(t(s_ind+1:e_ind),(X(1,s_ind+1:e_ind)/Ce0)*1e+3,'LineWidth',1.2)
ylabel('q/C_{E0} term (mV)')
xlabel('Time (seconds)')
title("Output voltage terms (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
set(gca,'FontSize',12)

subplot(3,1,3)
plot(t(s_ind+1:e_ind),X(5,s_ind+1:e_ind)*1e+3,'Color','b','LineWidth',1)
ylabel('Voltage (mV)')
xlabel('Time (seconds)')
title("Output Voltage (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
legend('State Space Model')
set(gca,'FontSize',12)

set(gcf,'Position',[100 50 900 750])

%% Plot spectra
figure
semilogx(fvec,20*log10(abs(X_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(f_ir,20*log10((abs(e_tf_freq(1:fsn/2))/fsn)*2),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+"dB re 20\muPa)")
legend('State Space Model','MATLAB IR Simulation')
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])

%%  Plot the output voltage only
s_ind = 1*fs;
e_ind = 1.5*fs;
maxssm = max(X(5,s_ind+1:e_ind)*1e+3);
figure
plot(t(s_ind+1:e_ind),X(5,s_ind+1:e_ind)*1e+3,'Color','b','LineWidth',1)
ylabel('Voltage (mV)')
xlabel('Time (seconds)')
title("Output Voltage (Input Signal = "...
+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+" dB re 20\muPa)")
legend('State Space Model')
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])

%% Transient analysis for different Rl

fs = 231500;    % Sampling frequency
A  = 1;           % Signal level
AdB = 20*log10(A/20e-6);  % Signal level in dB
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*fsig*t);
s_ind = 1*fs;
e_ind = T*fs;
Rl1 = 1000e+9;  % Load resistance

[X1,fvec1,X1_spec,THD1] = ssm_linear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,p_in,Rl1,s_ind,e_ind);

%% Plot spectra
Ohm = char(hex2dec('03A9'));
figure
semilogx(fvec,20*log10(abs(X_spec(5,:))),'Color','g','LineWidth',1)
hold on
semilogx(fvec,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(f_ir,20*log10((abs(e_tf_freq(1:fsn/2))/fsn)*2),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+"dB re 20\muPa)")
legend("State Space Model, R_L = 10G"+Ohm,"State Space Model, R_L = 1000G"+Ohm,'MATLAB IR Simulation')
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])