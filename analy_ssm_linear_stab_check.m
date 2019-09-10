% In this script the stability analysis for the linear state space model is
% carried out to determine the sampling frequency to be used

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

stabf = 1000:100:100e+6; % Frequencies to evaluate stability at
Rl = [1e+2 1e+4 1e+6 1e+8 1e+10 1e+12]; % Load resistances to evaluate stability at

% Call the function
eig_val = ssm_linear_stab_check(a,x0,E,Mmt,Cmt,Rmt,stabf,Rl);

%% Plot stability
Ohm = char(hex2dec('03A9'));
figure
subplot(1,2,1)
for ii = 1:3
    semilogx(stabf,log10(eig_val(ii,:)),'LineWidth',1.2)
    hold on
end
legend("R_L = 10^2"+Ohm,"R_L = 10^4"+Ohm,"R_L = 10^6"+Ohm)
xlabel('Sampling Frequency (Hz)')
ylabel('Maximum Eigen Value (log10)')
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
semilogx(stabf,log10(eig_val(4,:)),'LineWidth',1.2)
hold on
semilogx(stabf,log10(eig_val(5,:)),'--','LineWidth',1.2)
hold on
semilogx(stabf,log10(eig_val(6,:)),'-.','LineWidth',1.2)
legend("R_L = 10^8"+Ohm,"R_L = 10^{10}"+Ohm,"R_L = 10^{12}"+Ohm)
xlabel('Sampling Frequency (Hz)')
title('Stability Check','position',[12.4 3.08 0])
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 800 350])