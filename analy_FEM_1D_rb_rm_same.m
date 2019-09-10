% In this script, the static membrane deflection for different microphone
% parameters with the same membrane and backplate radius is calculated.

clc
clear
close all

rm = 10e-3;  % radius of the membrane
nn = 3001;  % number of nodes
E = 200;            % DC voltage
T = 3000;           % Tension of the membrane
d = 23.5*1e-6;      % distance between membrane and backplate
rb = rm;     % radius of the backplate
rw = 0;         % Width of the concentric ring
rp = 0;         % Position of the concentric ring from the center
func = "VA";    % Use Vamsy's function
e0 = 8.85*1e-12; % permittivity of air

%% FEM 1D Vs Analytical 1D when rb = rm

mP = "N";   % Apply a mean pressure
[w1,wuni1,ite1,wc1,rmp1,p1] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);
[w2,wuni2,ite2,wc2,rmp2] = analytical_1D(rm,nn,rb,E,d,T,mP);

mP = "Y";   % Apply a mean pressure
[w3,wuni3,ite3,wc3,rmp3,p3] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);
[w4,wuni4,ite4,wc4,rmp4] = analytical_1D(rm,nn,rb,E,d,T,mP);

%% Plot the membrane deflection

figure
% Uniform Pressure - Compare Analytical and FEM
subplot(2,2,1)
plot(rmp1*1e+3,wuni2*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(rmp2*1e+3,wuni1*1e+6,'--','LineWidth',1.2,'Color','b')
title('Analytical Vs FEM')
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
legend('Analytical Uniform','FEM Uniform')
grid minor
set(gca,'FontSize',12)

% Analytical Uni Vs FEM non Uni
Ca = e0*2*rmp3(2)*pi*sum(rmp3'./(d+wuni2));
Cf = e0*2*rmp1(2)*pi*sum(rmp1'./(d+w1));
subplot(2,2,2)
plot(rmp3*1e+3,wuni2*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(rmp1*1e+3,w1*1e+6,'LineWidth',1.2,'Color','g')
title('Analytical Vs FEM')
legend("Analytical Uniform, C_{E0}: "+num2str(Ca*1e+12,'%3.2f')+"pF",...
"FEM Nonuniform, C_{E0}: "+num2str(Cf*1e+12,'%3.2f')+"pF")
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
grid minor
set(gca,'FontSize',12)

% Plot the Non-uniform electrostatic pressure
subplot(2,2,3)
plot(rmp1*1e+3,p1,'LineWidth',1.2,'Color','r')
xlabel('Radius (mm)')
ylabel('Pressure (V^2F/m^3)')
legend('FEM Nonuniform')
title('Electrostatic Pressure')
grid minor
set(gca,'FontSize',12)

% Compare Non-averaged FEM and Averaged Analytical solutions
subplot(2,2,4)
plot(ite1,wc1*1e+6,'LineWidth',1.2,'Color','b')
xlabel('Number of Iterations')
ylabel('Deflection at r = 0 (\mum)')
title('Convergence Study')
legend('FEM Nonuniform')
grid minor
set(gca,'FontSize',12)

set(gcf,'position',[100 100 870 680]);

%% Number of Iteration Vs Center deflection
T  = fliplr(linspace(2000,3000,16));

for ii = 1:length(T)
    mP = "N";   % Apply a mean pressure
    [w,~,ite,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T(ii),mP,rw,rp,func);
    wc(ii) = w(1);
    nite(ii) = length(ite);
end

%% Plot the above results
figure
plot(abs(wc*1e+6),nite,'LineWidth',1.2,'Color','b')
ylabel('Number of Iterations')
xlabel('Deflection at r = 0 (\mum)')
title('Convergence Study')
legend('FEM Nonuniform')
grid minor
set(gca,'FontSize',12)
set(gcf,'position',[100 100 850 300]);
%% Compare FEM 1D and COMSOL
rb = rm;    % radius of the backplate
mP = "N";
[w1,wuni1,ite1,wc1,rmp1,p1] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);

% Load the COMSOL Model
model = mphload('2D_axisymm_st_rm_rb_same.mph');

Rmem = int2str(rm*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('Vpol',Vpol,'Polarization voltage');
model.param.set('Tm0',Tm0,'Membrane static tension');
model.param.set('Hm',Hm,'Air gap width');

% Disable the moving mesh components to get deflection due to uniform force
model.component('comp1').common('free1').active(false);
model.component('comp1').common('fix1').active(false);
model.component('comp1').common('disp1').active(false);
model.component('comp1').common('sym1').active(false);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2unidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',3);
w2uniEf = mpheval(model,"es.Ez",'edim','boundary','selection',3);

% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);
model.component('comp1').common('sym1').active(true);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2nunidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',3);
w2nuniEf = mpheval(model,"es.Ez",'edim','boundary','selection',3);
%% Compare FEM with COMSOL

% Find the maximum displacement from nonuniform def to set ylims
maxw = max(abs(w2nunidisp.d1));
maxy = -(maxw+0.5e-6)*1e+6;

figure
% Uniform force
subplot(2,2,1)
plot(rmp1*1e+3,wuni1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(w2unidisp.p(1,:)*1e+3,w2unidisp.d1*1e+6,'--','LineWidth',1.2,'Color','r')
ylabel('Deflection (\mum)')
xlabel('Radius (mm)')
lgd = legend('FEM Uniform','COMSOL Uniform');
lgd.FontSize = 10;
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Non-Uniform force
subplot(2,2,2)
plot(rmp1*1e+3,w1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(w2nunidisp.p(1,:)*1e+3,w2nunidisp.d1*1e+6,'--','LineWidth',1.2,'Color','r')
ylabel('Deflection (\mum)')
xlabel('Radius (mm)')
lgd = legend('FEM Nonuniform','COMSOL Nonuniform');
lgd.FontSize = 10;
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Uniform Electric field
subplot(2,2,3)
plot(w2uniEf.p(1,:)*1e+3,w2uniEf.d1,'LineWidth',1.2,'Color','b')
hold on
plot(w2nuniEf.p(1,:)*1e+3,w2nuniEf.d1,'--','LineWidth',1.2,'Color','r')
xlabel('Radius (mm)')
ylabel('Electric Field (V/m)')
title('COMSOL')
ylim([0.9*max(w2uniEf.d1) 1.1*max(w2nuniEf.d1)])
lgd = legend('Uniform','Nonuniform');
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)
set(gca,'position',[0.35 0.11 0.33 0.34])


% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Vs COMSOL 2D','t',supAxes);
set(gcf,'position',[50 50 870 680]);

%% Compare FEM with COMSOL - To plot extreme deflection cases

% Find the maximum displacement from nonuniform def to set ylims
maxw = max(abs(w1));
maxy = -(maxw+0.5e-6)*1e+6;

figure
% Uniform force
subplot(1,2,1)
plot(rmp1*1e+3,wuni1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(w2unidisp.p(1,:)*1e+3,w2unidisp.d1*1e+6,'--','LineWidth',1.2,'Color','r')
ylabel('Deflection (\mum)')
xlabel('Radius (mm)')
legend('FEM Uniform','COMSOL Uniform')
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Non-Uniform force
subplot(1,2,2)
plot(rmp1*1e+3,w1*1e+6,'LineWidth',1.2,'Color','b')
legend('FEM Nonuniform')
ylabel('Deflection (\mum)')
xlabel('Radius (mm)')
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.09 0.09 0.84 0.84];
suplabel('FEM 1D Vs COMSOL 2D','t',supAxes);
set(gcf,'position',[50 50 800 350]);