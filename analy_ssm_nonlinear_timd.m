% In this script, the transient analysis using a nonlinear state space model
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

%% Get the state variables and spectra - TIMD

fs = 231500;    % Sampling frequency
A  = 500;         % Signal level
A1dB = 20*log10(2*A/20e-6);
f1 = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*f1*t);
f2 = 1000;
fsig = 0;
p_in = p_in + A*sin(2*pi*f2*t);
mode = "TIMD";
Rl = 1e+10;
s_ind = 1*fs;
e_ind = T*fs;
qNL = "Y";
uNL = "Y";

[X1,fvec,X1_spec,TIMD1] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                        fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);
                                    
A = 2500;
A2dB = 20*log10(2*A/20e-6);
p_in = A*sin(2*pi*f1*t);
p_in = p_in + A*sin(2*pi*f2*t);

[X2,~,X2_spec,TIMD2] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                        fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);

%% Plot the membrane displacement and output voltage spectra
s_ind = 1*fs;
e_ind = 1.2*fs;
figure
subplot(2,1,1)
plot(t(s_ind+1:e_ind),X2(2,s_ind+1:e_ind)*1e+6,'Color','r','LineWidth',1.2)
ylabel('Displacement (\mum)')
xlabel('Time (seconds)')
title("Membrane Displacement (f1 = "+int2str(f1)+"Hz, f2 = "+int2str(f2)+"Hz)")
legend("Input: "+int2str(A2dB)+"dB")
set(gca,'FontSize',12)

% Find the DC shifts
maxp2 = max(X2(5,s_ind+1:e_ind));
maxn = max(-X2(5,s_ind+1:e_ind));
dcshift2 = maxp2-maxn;
subplot(2,1,2)
semilogx(fvec,20*log10(abs(X2_spec(5,:))),'Color','b','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (f1 = "+int2str(f1)+"Hz, f2 = "+int2str(f2)+"Hz)")
legend("Input: "+int2str(A2dB)+"dB, V_{peak}: "+num2str(maxp2,'%3.2f') +"V, DC shift: "+...
num2str(dcshift2,'%3.2f')+"V, TIMD: "+num2str(TIMD2(5),'%3.2f') + "%")
xlim([700 3000])
xticks([1000 2000])
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Get the state variables and spectra for diff f1 and f2

fs = 231500;    % Sampling frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
fsig = 0;
mode = "TIMD";
Rl = [10 50]*1e+9;
s_ind = 1*fs;
e_ind = T*fs;
qNL = "Y";
uNL = "Y";

A = [100 1000 5000 7000 9000 11000 13000 15000 16000 17000]/2;
for jj = 1:length(Rl)
    f1_1 = 20;
    f2_1 = 1000;
    for ii = 1:length(A)
        p_in = A(ii)*sin(2*pi*f1_1*t);
        p_in = p_in + A(ii)*sin(2*pi*f2_1*t);

        [~,~,~,TIMD1] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                                fs,fsig,f1_1,f2_1,p_in,mode,Rl(jj),s_ind,e_ind,qNL,uNL);
        timdA1(jj,ii) = TIMD1(5);
    end

    clear p_in
    f1_2 = 20;
    f2_2 = 5000;
    for ii = 1:length(A)
        p_in = A(ii)*sin(2*pi*f1_2*t);
        p_in = p_in + A(ii)*sin(2*pi*f2_2*t);

        [~,~,~,TIMD2] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                                fs,fsig,f1_2,f2_2,p_in,mode,Rl(jj),s_ind,e_ind,qNL,uNL);
        timdA2(jj,ii) = TIMD2(5);
    end


    clear p_in
    f1_3 = 100;
    f2_3 = 5000;
    for ii = 1:length(A)
        p_in = A(ii)*sin(2*pi*f1_3*t);
        p_in = p_in + A(ii)*sin(2*pi*f2_3*t);

        [~,~,~,TIMD3] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                                fs,fsig,f1_3,f2_3,p_in,mode,Rl(jj),s_ind,e_ind,qNL,uNL);
        timdA3(jj,ii) = TIMD3(5);
    end
end
%% Plot the above results
AdB = 20*log10(A/20e-6);
Ohm = char(hex2dec('03A9'));
figure
subplot(1,2,1)
plot(AdB,timdA1(1,:),'Color','r','LineWidth',1)
hold on
plot(AdB,timdA2(1,:),'--','Color','g','LineWidth',1)
hold on
plot(AdB,timdA3(1,:),'--','Color','b','LineWidth',1)
xlabel('Input Signal Level (dB re 20\muPa)')
ylabel('TIMD (%)')
title("Input Signal Level Vs TIMD (R_L = 10G"+Ohm+")")
legend("f1: "+int2str(f1_1)+"Hz, f2: "+int2str(f2_1)+"Hz","f1: "+int2str(f1_2)+"Hz, f2: "+int2str(f2_2)+"Hz",...
"f1: "+int2str(f1_3)+"Hz, f2: "+int2str(f2_3)+"Hz")
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
plot(AdB,timdA1(2,:),'Color','r','LineWidth',1)
hold on
plot(AdB,timdA2(2,:),'--','Color','g','LineWidth',1)
hold on
plot(AdB,timdA3(2,:),'--','Color','b','LineWidth',1)
xlabel('Input Signal Level (dB re 20\muPa)')
ylabel('TIMD (%)')
title("Input Signal Level Vs TIMD (R_L = 50G"+Ohm+")")
legend("f1: "+int2str(f1_1)+"Hz, f2: "+int2str(f2_1)+"Hz","f1: "+int2str(f1_2)+"Hz, f2: "+int2str(f2_2)+"Hz",...
"f1: "+int2str(f1_3)+"Hz, f2: "+int2str(f2_3)+"Hz")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 950 400])