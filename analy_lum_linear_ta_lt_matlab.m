% In this script the transient responses from LTSpice and MATLAB
% simulation are compared. Remember to not consider the transient
% response in LTSpice before viewing the FFT.

clc
clear
close all

% Load the output from LTspice
fsig = 1000;
Pi = 1000;

% Load the transient response saved from LTSpice
% filename = "4133_Linear_TA_" + int2str(fsig) + "Hz_"+int2str(Pi) + "V_TA.txt";
% [time,amp]=read_spice_time(filename);

% Load the frequency response saved from LTSpice
filename = "4133_Linear_TA_" + int2str(fsig) + "Hz_"+int2str(Pi) + "V_FFT.txt";
[freq,mag,phase]=read_spice_freq_reim(filename);

% Calculate THD from the FFT
Y = 10.^(mag/20);
nmax = 20; % Maximum number of harmonics to consider for THD calculation
Nsl = 1;      % Maximum number of side libes to consider for THD calculation
[THD_lt,~] = thd_nu(Y,fsig,freq,nmax,Nsl);

% Get the frequency response with the same number of points of that LTSpice
% output
P_in = 1;
[~,e_noref,~,~,~,~] = mic_tf_sim(P_in,freq);
PdB = 20*log10(Pi/20e-6);

% Apply ifft to get Impulse Response
E_noref = [e_noref' fliplr(e_noref')];
mic_ir = ifft(E_noref,'symmetric'); 
fs = 2*ceil(max(freq));
tr = 0:1/fs:2-1/fs;
f_ir = 0:fs/2 - 1;

% Convolution of IR with input signal
inp = Pi*sin(2*pi*fsig*tr);
e_tf_amp = conv(inp,mic_ir);     % Transient response of the TF output
e_tf_freq = fft(e_tf_amp(fs+1:2*fs));  % Only take till the length of input signal
Y = abs(e_tf_freq(1:fs/2));
[THD_tf,~] = thd_nu(Y,fsig,f_ir,nmax,Nsl);

% Get the time scale for the MATLAB simulated output
T = length(e_tf_amp)/fs;
tn = 0:1/fs:T-1/fs;

%% Plot the time domain and FFT response

% Find the start and end indices to plot for LTSpice
indS = find(time >= 1,1);
indE = find(time >= 1.1,1);
figure
subplot(2,1,1)
plot(time(indS:indE),amp(indS:indE)*1e+3,'LineWidth',1,'Color','b')
hold on
% Find the start and end indices to plot for MATLAB output
indS = find(tn >= 1,1);
indE = find(tn >= 1.1,1);
plot(tn(indS:indE),e_tf_amp(indS:indE)*1e+3,'--','LineWidth',1,'Color','r')
xlabel('Time (seconds)')
ylabel('Voltage (mV)')
title("Output Voltage (Input Signal = "+int2str(fsig)+" Hz puretone at "+int2str(PdB)+"dB re 20\muPa)")
legend('LTspice Simulation','MATLAB IR Simulation')
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(freq,mag,'LineWidth',1,'Color','b')
hold on
semilogx(f_ir,20*log10((abs(e_tf_freq(1:fs/2))/fs)*2),'--','LineWidth',1,'Color','r')
grid minor
xlabel('Frequency (Hz)')
ylabel('Sensitivity (dB re 1V/Pa)')
title("Frequency Spectrum (Input Signal = "+int2str(fsig)+" Hz puretone at "+int2str(PdB)+"dB re 20\muPa)")
legend("LTspice Simulation, THD: " + num2str(THD_lt,'%3.2f')+"%","MATLAB IR Simulation,THD:"+ num2str(THD_tf,'%3.2f')+"%")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Only plot the spectrum for high frequencies as the Spice has very sampling rate
figure
semilogx(freq,mag,'LineWidth',1,'Color','b')
hold on
semilogx(f_ir,20*log10((abs(e_tf_freq(1:fs/2))/fs)*2),'--','LineWidth',1,'Color','r')
grid minor
xlabel('Frequency (Hz)')
ylabel('Sensitivity (dB re 1V/Pa)')
title("Frequency Spectrum (Input Signal = "+int2str(fsig)+" Hz puretone at "+int2str(PdB)+"dB re 20\muPa)")
legend("LTspice Simulation, THD: " + num2str(THD_lt,'%3.2f')+"%","MATLAB IR Simulation,THD:"+ num2str(THD_tf,'%3.2f')+"%")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])