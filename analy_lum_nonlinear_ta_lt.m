% In this script the transient responses from LTSpice and MATLAB
% simulation are compared. Remember to not consider the transient
% response in LTSpice before viewing the FFT.

clc
clear
close all

% Load the LTspice output
fsig = 20;
[time1,amp1]=read_spice_time("4133_Nonlinear_Bv_TA_20Hz_1V_TA.txt");
[time2,amp2]=read_spice_time("4133_Nonlinear_Bv_TA_20Hz_100V_TA.txt");

% Get the frequency response with the same number of points of that LTSpice
% output
freq = linspace(0,50000,50001);
P_in = 1;
[~,e_noref,~,~,~,~] = mic_tf_sim(P_in,freq);

% Apply ifft to get Impulse Response
E_noref = [e_noref fliplr(e_noref)];
mic_ir = ifft(E_noref,'symmetric'); 
fs = 2*ceil(max(freq));
tr = 0:1/fs:2-1/fs;

% Convolution of IR with input signal
inp = P_in*sin(2*pi*fsig*tr);
e_tf_amp = conv(inp,mic_ir);     % Transient response of the TF output

% Get the time scale for the MATLAB simulated output
T = length(e_tf_amp)/fs;
tn = 0:1/fs:T-1/fs;

%% Plot and compare
PdB = 94;
% Find the start and end indices to plot for LTSpice
indS = find(time1 >= 1,1);
indE = find(time1 >= 1.1,1);
figure
subplot(2,1,1)
plot(time1(indS:indE),amp1(indS:indE)*1e+3,'LineWidth',1,'Color','b')
hold on
% Find the start and end indices to plot for MATLAB output
indS = find(tn >= 1,1);
indE = find(tn >= 1.1,1);
plot(tn(indS:indE),e_tf_amp(indS:indE)*1e+3,'--','LineWidth',1,'Color','r')
xlabel('Time (seconds)')
ylabel('Voltage (mV)')
legend('LTspice Simulation','MATLAB IR Simulation')
title("Transient response (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(PdB)+"dB re 20\muPa)")
set(gca,'FontSize',12)

PdB = 134;
subplot(2,1,2)
plot(time2,amp2*1e+3,'LineWidth',1,'Color','b')
ylabel('Voltage (mV)')
xlabel('Time (seconds)')
title("Transient response (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(PdB)+"dB re 20\muPa)")
set(gca,'FontSize',12)
legend('LTspice Simulation')
set(gcf,'Position',[100 100 900 500])

%%  Caclculate THD at high levels only including BMechNL

clc
clear
close all

% Load the LTspice output
fsig = 1000;
[time1,amp1]=read_spice_time("4133_Nonlinear_Bv_MechNL_TA_20Hz_100V_TA.txt");
[time2,amp2]=read_spice_time("4133_Nonlinear_Bv_MechNL_TA_20Hz_1000V_TA.txt");

[freq1,mag1,phase1]=read_spice_freq_reim("4133_Nonlinear_Bv_MechNL_TA_1000Hz_100V_FFT.txt");
[freq2,mag2,phase2]=read_spice_freq_reim("4133_Nonlinear_Bv_MechNL_TA_1000Hz_1000V_FFT.txt");

% Calculate THD from the FFT for 134dB
Y = 10.^(mag1/20);
nmax = 10;   % Maximum number of harmonics to consider for THD calculation
Nsl = 1;      % Maximum number of side libes to consider for THD calculation
[THD1,~] = thd_nu(Y,fsig,freq1,nmax,Nsl);

% Calculate THD from the FFT for 154dB
Y = 10.^(mag2/20);
nmax = 10;   % Maximum number of harmonics to consider for THD calculation
Nsl = 1;      % Maximum number of side libes to consider for THD calculation
[THD2,~] = thd_nu(Y,fsig,freq2,nmax,Nsl);

%% Plot and compare
% Find the start and end indices to plot for LTSpice
indS = find(time1 >= 1,1);
indE = find(time1 >= 1.1,1);
figure
subplot(2,1,1)
plot(time1(indS:indE),amp1(indS:indE),'LineWidth',1,'Color','b')
hold on
% Find the start and end indices to plot for LTSpice
indS = find(time2 >= 1,1);
indE = find(time2 >= 1.1,1);
plot(time2(indS:indE),amp2(indS:indE),'--','LineWidth',1,'Color','r')
xlabel('Time (seconds)')
ylabel('Voltage (V)')
legend('LTspice Simulation, Input: 134dB','LTspice Simulation, Input: 154dB')
title("Transient response (Input Signal = "+int2str(fsig)+"Hz Puretone)")
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(freq1,mag1,'LineWidth',1,'Color','b')
hold on
semilogx(freq2,mag2,'--','LineWidth',1,'Color','r')
grid minor
xlabel('Frequency (Hz)')
ylabel('Sensitivity (dB re 1V/Pa)')
title("Frequency Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("LTspice Simulation, Input: 134dB, THD: "+num2str(THD1,'%3.2f')+"%","LTspice Simulation, Input: 154dB, THD: "+num2str(THD2,'%3.2f')+"%")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Only plot the spectrum for higher frequencies

figure
semilogx(freq1,mag1,'LineWidth',1,'Color','b')
hold on
semilogx(freq2,mag2,'--','LineWidth',1,'Color','r')
grid minor
xlabel('Frequency (Hz)')
ylabel('Sensitivity (dB re 1V/Pa)')
title("Frequency Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("LTspice Simulation, Input: 134dB, THD: "+num2str(THD1,'%3.2f')+"%","LTspice Simulation, Input: 154dB, THD: "+num2str(THD2,'%3.2f')+"%")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])
