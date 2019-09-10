clc
clear
close all

rm = 10e-3;  % radius of the membrane
nn = 3001;  % number of nodes
mode = 'r'; % r to compute only across the radius and d to compute across diameter
E = 200;            % DC voltage
T = 3000;           % Tension of the membrane
d = 23.5*1e-6;      % distance between membrane and backplate
d1 = 2e-3;      % slit gap
rb = rm - d1;    % radius of the backplate
rw = 0;         % Width of the concentric ring
rp = 0;         % Position of the concentric ring from the center
func = "VA";    % Use Vamsy's function
tm = 7e-6;      % Membrane thickness
rhom = 8300;    % Membrane density
surfd = rhom*tm; % Membrane surface density
e0 = 8.85*1e-12; % permittivity of air
Rl = 100e+9;

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
w3 = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6]);
fC = model.study('std1').feature('frlin').getDoubleArray('plist');

%% FEM 1D

mP = "N";   % Apply a mean pressure
f = [fC(1) fC(6) fC(11) fC(16)];         % Input signal frequency
pdb = 20*log10(pin/20e-6);

% Get the static deflection first
[w1,wuni1,ite1,wc1,rmp1,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);
ind = find(rmp1 <= rb);
Cs = e0*2*pi*rmp1(2)*sum(rmp1(ind)'./(d+w1(ind)));
Rl = 1e+9;
RCt = Cs*Rl;

for ii = 1:4

    % Use the static as initial condition and calculate the dynamic deflection
    [wt,rmp2,~] = FEM_1D_dyn_cyl(rm,nn,rb,E,d,T,func,f(ii),pin,surfd,w1);
    w2(ii,:) = wt(2,:);
end

%% Plot the membrane deflection

figure
subplot(2,2,1)
plot(rmp2*1e+3,w2(1,:)*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(1,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(1))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,2)
plot(rmp2*1e+3,w2(2,:)*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(6,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(2))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,3)
plot(rmp2*1e+3,w2(3,:)*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(11,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(3))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,4)
yyaxis right
plot(rmp2*1e+3,w2(4,:)*1e+6,'LineWidth',1.2)
ylabel('Dynamic deflection (\mum)')
yyaxis left
plot(w3.p(1,:)*1e+3,abs(w3.d1(16,:))*1e+6,'--','LineWidth',1.2)
title("f = "+int2str(f(4))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('COMSOL 2D','FEM 1D')
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .09 0.84 0.84];
%suplabel('Membrane radial direction (mm)','x',supAxes);
%suplabel('Dyanmic deflection only (\mum)','y',supAxes);
%suplabel('FEM 1D Vs COMSOL 2D','t',supAxes);
set(gcf,'position',[50 50 870 680]);

%% Plot the membrane deflection

figure
subplot(2,2,1)
plot(rmp2*1e+3,-(abs(w2(1,:))-abs(w1'))*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(1,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(1))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,2)
plot(rmp2*1e+3,-(abs(w2(2,:))-abs(w1'))*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(6,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(2))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,3)
plot(rmp2*1e+3,-(abs(w2(2,:))-abs(w1'))*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(w3.p(1,:)*1e+3,-abs(w3.d1(11,:))*1e+6,'--','LineWidth',1.2,'Color','b')
title("f = "+int2str(f(3))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('FEM 1D','COMSOL 2D')
grid minor
set(gca,'FontSize',12)

subplot(2,2,4)
yyaxis right
plot(rmp2*1e+3,-(abs(w2(4,:))-abs(w1'))*1e+6,'LineWidth',1.2)
ylabel('Dynamic deflection (\mum)')
yyaxis left
plot(w3.p(1,:)*1e+3,abs(w3.d1(16,:))*1e+6,'--','LineWidth',1.2)
title("f = "+int2str(f(4))+"Hz")
xlabel('Radius (mm)')
ylabel('Dynamic deflection (\mum)')
legend('COMSOL 2D','FEM 1D')
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .09 0.84 0.84];
%suplabel('Membrane radial direction (mm)','x',supAxes);
%suplabel('Dyanmic deflection only (\mum)','y',supAxes);
%suplabel('FEM 1D Vs COMSOL 2D','t',supAxes);
set(gcf,'position',[50 50 870 680]);