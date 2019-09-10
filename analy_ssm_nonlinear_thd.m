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

%% Get the state variables and validate against Linear model

fs = 2315000;    % Sampling frequency
A  = 1;         % Signal level
AdB = 20*log10(A/20e-6);
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*fsig*t);
mode = "THD";
f1 = 0;
f2 = 0;
s_ind = 1*fs;
e_ind = T*fs;
Rl = 10e+9;
qNL = "Y";
uNL = "Y";

% Nonlinear model
tic
[X1,fvec1,X1_spec,THD1] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);
etaNL = toc;

% Linear model
tic
[X2,fvec2,X2_spec,THD2] = ssm_linear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,p_in,Rl,s_ind,e_ind);
etaL = toc;
%% Plot and validate the output voltage time response and spectra
s_ind = 1*fs;
e_ind = 1.5*fs;
figure
subplot(2,1,1)
plot(t(s_ind+1:e_ind),X2(5,s_ind+1:e_ind)*1e+3,'Color','b','LineWidth',1.2)
hold on
plot(t(s_ind+1:e_ind),X1(5,s_ind+1:e_ind)*1e+3,'--','Color','r','LineWidth',1.2)
ylabel('Voltage (mV)')
xlabel('Time (seconds)')
title("Output Voltage (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+"dB re 20\muPa)")
legend('Linear State Space Model','Nonlinear State Space Model')
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(fvec2,20*log10(abs(X2_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec1,20*log10(abs(X1_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(AdB)+"dB re 20\muPa)")
legend('Linear State Space Model','Nonlinear State Space Model')
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Get the state variables

fs = 2315000;    % Sampling frequency
A  = 1;         % Signal level
A1dB = 20*log10(A/20e-6);
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
p_in = A*sin(2*pi*fsig*t);
mode = "THD";
f1 = 0;
f2 = 0;
s_ind = 1*fs;
e_ind = T*fs;
Rl = 50e+9;
qNL = "N";
uNL = "N";

% Nonlinear model for A1
[X1,fvec1,X1_spec,THD1] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);

A  = 100;         % Signal level
A2dB = 20*log10(A/20e-6);
p_in = A*sin(2*pi*fsig*t);
% Nonlinear model for A2
[X2,fvec2,X2_spec,THD2] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);

%% Plot and validate the output voltage time response and spectra
s_ind = 1*fs;
e_ind = 1.1*fs;
maxp1 = max(X1(5,s_ind+1:e_ind));
maxn = max(-X1(5,s_ind+1:e_ind));
dcshift1 = maxp1-maxn;
maxp2 = max(X2(5,s_ind+1:e_ind));
maxn = max(-X2(5,s_ind+1:e_ind));
dcshift2 = maxp2-maxn;
figure
subplot(2,1,1)
plot(t(s_ind+1:e_ind),X1(5,s_ind+1:e_ind),'Color','b','LineWidth',1.2)
hold on
plot(t(s_ind+1:e_ind),X2(5,s_ind+1:e_ind),'--','Color','r','LineWidth',1.2)
ylabel('Voltage (V)')
xlabel('Time (seconds)')
title("Output Voltage (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Nonlinear State Space Model, Input: "+int2str(A1dB)+"dB, DC Shift: "+num2str(dcshift1,'%3.2f')+"V",...
"Nonlinear State Space Model, Input: "+int2str(A2dB)+"dB, DC Shift: "+num2str(dcshift2,'%3.2f')+"V")
set(gca,'FontSize',12)

subplot(2,1,2)
semilogx(fvec2,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec1,20*log10(abs(X2_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Nonlinear State Space Model, Input: "+int2str(A1dB)+"dB, THD: "+num2str(THD1(5),'%3.2f') + "%",...
    "Nonlinear State Space Model, Input: "+int2str(A2dB)+"dB, THD: "+num2str(THD2(5),'%3.2f') + "%")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Plot only spectra
figure
semilogx(fvec2,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec1,20*log10(abs(X2_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Nonlinear State Space Model, Input: "+int2str(A1dB)+"dB, THD: "+num2str(THD1(5),'%3.2f') + "%",...
    "Nonlinear State Space Model, Input: "+int2str(A2dB)+"dB, THD: "+num2str(THD2(5),'%3.2f') + "%")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])

%% Plot the membrane displacement and output voltage spectra
s_ind = 1*fs;
e_ind = 1.1*fs;
figure
subplot(2,1,1)
plot(t(s_ind+1:e_ind),X1(2,s_ind+1:e_ind)*1e+6,'Color','b','LineWidth',1.2)
hold on
plot(t(s_ind+1:e_ind),X2(2,s_ind+1:e_ind)*1e+6,'--','Color','r','LineWidth',1.2)
ylabel('Displacement (\mum)')
xlabel('Time (seconds)')
title("Membrane Displacement (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Input: "+int2str(A1dB)+"dB","Input: "+int2str(A2dB)+"dB")
set(gca,'FontSize',12)

% Find the DC shifts
maxp1 = max(X1(5,s_ind+1:e_ind));
maxn = max(-X1(5,s_ind+1:e_ind));
dcshift1 = maxp1-maxn;
maxp2 = max(X2(5,s_ind+1:e_ind));
maxn = max(-X2(5,s_ind+1:e_ind));
dcshift2 = maxp2-maxn;
subplot(2,1,2)
semilogx(fvec1,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec2,20*log10(abs(X2_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Input: "+int2str(A1dB)+"dB, V_{peak}: "+num2str(maxp1,'%3.2f') +"V, DC shift: "+num2str(dcshift1,'%3.2f')+"V, THD: "+num2str(THD1(5),'%3.2f') + "%",...
    "Input: "+int2str(A2dB)+"dB, V_{peak}: "+num2str(maxp2,'%3.2f') +"V, DC shift: "+num2str(dcshift2,'%3.2f')+"V, THD: "+num2str(THD2(5),'%3.2f') + "%")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 500])

%% Plot only the output voltage spectra
% Find the DC shifts
maxp1 = max(X1(5,s_ind+1:e_ind));
maxn = max(-X1(5,s_ind+1:e_ind));
dcshift1 = maxp1-maxn;
maxp2 = max(X2(5,s_ind+1:e_ind));
maxn = max(-X2(5,s_ind+1:e_ind));
dcshift2 = maxp2-maxn;
figure
semilogx(fvec2,20*log10(abs(X1_spec(5,:))),'Color','b','LineWidth',1)
hold on
semilogx(fvec1,20*log10(abs(X2_spec(5,:))),'--','Color','r','LineWidth',1)
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title("Output Voltage Spectrum (Input Signal = "+int2str(fsig)+"Hz Puretone)")
legend("Input: "+int2str(A1dB)+"dB, V_{peak}: "+num2str(maxp1,'%3.2f') +"V, DC shift: "+num2str(dcshift1,'%3.2f')+"V, THD: "+num2str(THD1(5),'%3.2f') + "%",...
    "Input: "+int2str(A2dB)+"dB, V_{peak}: "+num2str(maxp2,'%3.2f') +"V, DC shift: "+num2str(dcshift2,'%3.2f')+"V, THD: "+num2str(THD2(5),'%3.2f') + "%")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])

%% Compute output for multiple input levels - THD Vs Membrane displacement

fs = 2315000;    % Sampling frequency
fsig = 20;      % Signal frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
mode = "THD";
f1 = 0;
f2 = 0;
s_ind = 1*fs;
e_ind = T*fs;
Rl = 10e+9;
qNL = "Y";
uNL = "Y";

% Run in a loop
A = [100 1000 5000 7000 9000 11000 13000 15000 16000 17000];
for ii = 1:length(A)
    p_in = A(ii)*sin(2*pi*fsig*t);
    [X,fvec,X_spec,THD] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,f1,f2,p_in,mode,Rl,s_ind,e_ind,qNL,uNL);
    
    % Save the peak membrane displacement
    mDisp(ii) = max(X(2,s_ind+1:e_ind));
    
    % Save the THD
    thdA(ii) = THD(5);
end

% Run in a loop
AR = 16000;
p_in = AR*sin(2*pi*fsig*t);
Rl = [10 20 30 40 50]*1e+9;
for jj = 1:length(Rl)
    
    [X,fvec,X_spec,THD] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                       fs,fsig,f1,f2,p_in,mode,Rl(jj),s_ind,e_ind,qNL,uNL);
    
    % Save the THD
    thdRl(jj) = THD(5);
end

%% Plot the displacement and THD from above
AdB = 20*log10(abs(A)/20e-6);
ARdB = 20*log10(abs(AR)/20e-6);
Ohm = char(hex2dec('03A9'));
figure
subplot(1,2,1)
yyaxis left
plot(AdB,mDisp*1e+6,'--','LineWidth',1)
ylabel('Membrane Displacement (\mum)')
yyaxis right
plot(AdB,thdA,'*')
ylabel('THD (%)')
xlabel('Input Signal Level (dB re 20\muPa)')
title("Input Signal = "+int2str(fsig)+"Hz Puretone")
grid minor
set(gca,'FontSize',12)
subplot(1,2,2)
plot(Rl*1e-9,thdRl,'*','Color','r')
xlabel("Load Resistance (G"+Ohm+")")
ylabel('THD (%)')
title("Input Signal = "+int2str(fsig)+"Hz Puretone at "+int2str(ARdB)+"dB re 20\muPa")
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 950 400])

%% Compute output for multiple frequencies

fs = 3000000;    % Sampling frequency
T = 2;          % length of signal in seconds
t = 0:1/fs:T-1/fs;
mode = "THD";
f1 = 0;
f2 = 0;
s_ind = 1*fs;
e_ind = T*fs;
qNL = "Y";
uNL = "Y";

% Run in a loop
A = 16000;
Rl = [10 50]*1e+9;
fsig = [10 20 100 1000 5000 8000 12000 16000 20000];
for ii = 1:length(Rl)
    for jj = 1:length(fsig)
    
        p_in = A*sin(2*pi*fsig(jj)*t);
        [~,~,~,THD] = ssm_nonlinear(a,x0,E,Mmt,Cmt,Rmt,...
                                           fs,fsig(jj),f1,f2,p_in,mode,Rl(ii),s_ind,e_ind,qNL,uNL);

        % Save the THD
        thdfreq(ii,jj) = THD(5);
    end
end

%% Plot the above results
Ohm = char(hex2dec('03A9'));
AdB = 20*log10(A/20e-6);
figure
semilogx(fsig,thdfreq(1,:),'Color','b','LineWidth',1)
hold on
semilogx(fsig,thdfreq(2,:),'--','Color','r','LineWidth',1)
ylabel('THD (%)')
xlabel('Frequency (Hz)')
title("Frequency Vs THD (Input Signal Level = "+int2str(AdB)+"dB re 20\muPa)")
legend("Load Resistance: "+int2str(Rl(1)*1e-9)+ "G"+Ohm,...
    "Load Resistance: "+int2str(Rl(2)*1e-9)+ "G"+Ohm)
grid minor
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 900 250])