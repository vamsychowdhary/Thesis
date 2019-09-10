% In this script the linear transfer function from LTSpice and MATLAB
% simulation is compared.

clc
clear
close all

% Get the Measured and calibrated response
[mic_4133,mic_4133_ph,f_4133] = getMeasResp4133;

% Load the LTspice output
[f_lt,mag_lt,ph_lt]=read_spice_freq_reim("4133_Linear_AC.txt");

% Get the TF simulated in MATLAB
Pi = 1;
[e_ref,~,ph_ref,~,PdB] = mic_tf_sim(Pi,f_lt);
e_ref_dB = 20*log10(abs(e_ref));

% Unwrap LTSpice phase
up_lt = unwrap(ph_lt);
up_lt = up_lt*180/pi;
up_lt = up_lt + 360;

% Plot the transfer functions including diaphragm reflections
figure
subplot(2,1,1)
semilogx(f_lt,mag_lt,'LineWidth',1.5,'Color','b')
hold on
semilogx(f_lt,e_ref_dB,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133,'--','LineWidth',1.5,'Color','g')
grid minor
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Frequency Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('LTspice Simulation','MATLAB TF Simulation','Measured Response')
ylim([-55 -36])
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(f_lt,up_lt,'LineWidth',1.5,'Color','b')
hold on
semilogx(f_lt,ph_ref,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133_ph + 180,'--','LineWidth',1.5,'Color','g')
grid minor
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title("Phase Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('LTspice Simulation','MATLAB TF Simulation','Measured Response')
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Plot the frequency response by excluding the diaphragm reflections from
% TF simulation

% Get the Measured and calibrated response
[mic_4133,mic_4133_ph,f_4133] = getMeasResp4133;

% Get the TF simulated in MATLAB
Pi = 1;
[e_ref,e_noref,ph_ref,ph_noref,PdB,Ts] = mic_tf_sim(Pi,f_lt);
e_ref_dB = 20*log10(abs(e_ref));
e_noref_dB = 20*log10(abs(e_noref));

% Plot the transfer functions including diaphragm reflections
figure
subplot(2,1,1)
semilogx(f_lt,e_ref_dB,'LineWidth',1.5,'Color','b')
hold on
semilogx(f_lt,e_noref_dB,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133,'--','LineWidth',1.5,'Color','g')
grid minor
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Frequency Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('TF Simulation including reflections','TF Simulation excluding reflections',...
'Measured Response')
ylim([-60 -36])
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(f_lt,ph_ref,'LineWidth',1.5,'Color','b')
hold on
semilogx(f_lt,ph_noref,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133_ph + 180,'--','LineWidth',1.5,'Color','g')
grid minor
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title("Phase Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('TF Simulation including reflections','TF Simulation excluding reflections',...
'Measured Response')
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

% Plot the diaphragm reflections
figure
semilogx(f_lt,20*log10(abs(Ts)),'LineWidth',1.5,'Color','b')
grid minor
ylabel('Gain (dB)')
xlabel('Frequency (Hz)')
title("Membrane Reflections Transfer Function")
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])

%% Plot the frequency response by excluding the diaphragm reflections from
% TF simulation and tweaking Ras and Mas to match the simulated response
% with measured. Please change inside mic_tf_sim function manually.

% Get the Measured and calibrated response
[mic_4133,mic_4133_ph,f_4133] = getMeasResp4133;

% Get the TF simulated in MATLAB
Pi = 1;
[~,e_noref,~,ph_noref,PdB,Ts] = mic_tf_sim(Pi,f_lt);
e_noref_dB = 20*log10(abs(e_noref));

% Plot the transfer functions including diaphragm reflections
figure
subplot(2,1,1)
semilogx(f_lt,e_noref_dB,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133,'--','LineWidth',1.5,'Color','g')
grid minor
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Frequency Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('TF Simulation excluding reflections','Measured Response')
ylim([-60 -36])
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(f_lt,ph_noref,'--','LineWidth',1.5,'Color','r')
hold on
semilogx(f_4133,mic_4133_ph + 180,'--','LineWidth',1.5,'Color','g')
grid minor
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title("Phase Response (Input Pressure = "+int2str(PdB)+"dB re 20\muPa)")
legend('TF Simulation excluding reflections','Measured Response')
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])