clc
clear
close all

rm = 10e-3;  % radius of the membrane
nn = 3001;  % number of nodes
mode = 'r'; % r to compute only across the radius and d to compute across diameter
E = 200;            % DC voltage
T = 3000;           % Tension of the membrane
d = 100*1e-6;      % distance between membrane and backplate
d1 = 2e-3;      % slit gap
rb = rm - d1;    % radius of the backplate
rw = 0;         % Width of the concentric ring
rp = 0;         % Position of the concentric ring from the center
func = "VA";    % Use Vamsy's function
tm = 7e-6;      % Membrane thickness
rhom = 8300;    % Membrane density
surfd = rhom*tm; % Membrane surface density
e0 = 8.85*1e-12; % permittivity of air
Rl = 1000e+9;

%% Load COMSOL model
pin = 1;        % Input signal pressure level

% Load the COMSOL Model
model = mphload('2D_axisymm_dy_rm_rb_diff1.mph');

Pin = int2str(pin)+ "[Pa]";
Hm = int2str(d*1e+6)+ "[um]";
G = num2str(d1*1e+3,'%3.2f')+ "[mm]";
Rmem = int2str(rm*1e+3)+ "[mm]";
Vpol = int2str(E)+ "[V]";
Rpreamp = int2str(Rl*1e-9)+ "[Gohm]";
Tm0 = int2str(T)+ "[N/m]";

% Update model parameters accordingly
model.param.set('pin', Pin, 'External incident pressure');
model.param.set('Rmem', Rmem, 'Membrane Radius');
model.param.set('Hm', Hm, 'Air gap thickness');
model.param.set('Vpol', Vpol, 'Polarization voltage');
model.param.set('G', G , 'Slit Gap width');
model.param.set('Rpreamp', Rpreamp, 'Preamplifier input impedance');
model.param.set('Tm0', Tm0, 'Tension');

% Run the static part
model.study("std1").run

% Load the displacement z component and e field to a variable
eCdB = mpheval(model,"es.V0_1");
fC = model.study('std1').feature('frlin').getDoubleArray('plist');

%% FEM 1D
mP = "N";   % Apply a mean pressure
pdb = 20*log10(pin/20e-6);

% Get the static deflection first
[w1,wuni1,ite1,wc1,rmp1,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);

ind = find(rmp1 <= rb);
Cs = e0*2*pi*rmp1(2)*sum(rmp1(ind)'./(d+w1(ind)));
RCt = Cs*Rl;

for ii = 1:length(fC)
    
    % Use the static as initial condition and calculate the dynamic deflection
    [~,~,~,e,edB] = FEM_1D_dyn_cyl(rm,nn,rb,E,d,T,func,fC(ii),pin,surfd,w1);
    ef(1,ii) = e(1);
    ef(2,ii) = e(2);
    edBf(1,ii) = edB(1);
    edBf(2,ii) = edB(2);
end 

%% Plot the frequency response

% Unwrap the phase
up = unwrap(atan2(imag(eCdB.d1(:,1)),real(eCdB.d1(:,1))));
%up = atan2(imag(eCdB.d1(:,1)),real(eCdB.d1(:,1)));
%up = angle(eCdB.d1(:,1)');
up = up*180/pi;
%up = up + 360;

figure
subplot(2,1,1)
semilogx(fC,edBf(1,:),'.-','LineWidth',1.2,'Color','r')
hold on
semilogx(fC,edBf(2,:),'LineWidth',1.2,'Color','g')
hold on
semilogx(fC,20*log10(abs(eCdB.d1(:,1))),'--','LineWidth',1.2,'Color','b')
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title('Frequency Response')
legend('FEM 1D Static','FEM 1D No Static','COMSOL 2D Static')
grid minor
set(gca,'FontSize',12)
% set(gcf,'position',[100 50 900 270]);

subplot(2,1,2)
semilogx(fC,up,'LineWidth',1.2,'Color','b')
ylabel('Phase (degress)')
xlabel('Frequency (Hz)')
title('Unwrapped Phase Response (COMSOL 2D)')
grid minor
set(gca,'FontSize',12)
set(gcf,'position',[100 50 900 550]);

%%
figure
% subplot(2,1,1)
semilogx(fC,edBf(1,:),'LineWidth',1.2,'Color','r')
hold on
semilogx(fC,20*log10(abs(eCdB.d1(:,1))),'--','LineWidth',1.2,'Color','b')
ylabel('Sensitivity (dB re 1V/Pa)')
xlabel('Frequency (Hz)')
title('Frequency Response')
legend('FEM 1D Static','COMSOL 2D Static')
grid minor
set(gca,'FontSize',12)
set(gcf,'position',[100 50 900 270]);