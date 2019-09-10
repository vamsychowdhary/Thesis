clc
clear
close all

rm = 10e-3;  % radius of the membrane
nn = 3001;  % number of nodes
mode = 'r'; % r to compute only across the radius and d to compute across diameter
E = 200;            % DC voltage
T = 3000;           % Tension of the membrane
d = 23.5*1e-6;      % distance between membrane and backplate
rb = rm - 2e-3;    % radius of the backplate
rw = 0;         % Width of the concentric ring
rp = 0;         % Position of the concentric ring from the center
func = "VA";    % Use Vamsy's function

%% FEM 1D Vs Analytical 1D when rb ~= rm

mP = "N";   % Apply a mean pressure
[w1,wuni1,ite1,wc1,rmp1,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);

mP = "Y";   % Apply a mean pressure
[w2,wuni2,ite2,wc2,rmp2,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);
[w3,wuni3,ite3,wc3,rmp3] = analytical_1D(rm,nn,rb,E,d,T,mP);

%% Plot the membrane deflection

figure
% Uniform force
subplot(2,2,1)  
plot(rmp1*1e+3,wuni1*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(rmp3*1e+3,wuni3*1e+6,'--','LineWidth',1.2,'Color','b')
title('Uniform force')
legend('FEM','Analytical')
grid minor
set(gca,'FontSize',12)

% Compare Averaged FEM and Analytical solutions
subplot(2,2,2)
plot(rmp2*1e+3,w2*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(rmp3*1e+3,w3*1e+6,'--','LineWidth',1.2,'Color','b')
title('Non-uniform force')
legend('Averaged FEM','Averaged Analytical')
grid minor
set(gca,'FontSize',12)

% Compare Non-averaged FEM and Aevraged Analytical solutions
subplot(2,2,3)
plot(rmp1*1e+3,w1*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(rmp3*1e+3,w3*1e+6,'--','LineWidth',1.2,'Color','b')
title('Non-uniform force')
legend('Non-averaged FEM','Averaged Analytical')
grid minor
set(gca,'FontSize',12)
set(gca,'position',[0.35 0.11 0.33 0.34])

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .09 0.8 0.84];
suplabel('Membrane radial direction (mm)','x',supAxes);
suplabel('Membrane deflection (\mum)','y',supAxes);
set(gcf,'position',[100 150 800 600]);

%% Plot the convergence between FEM 1D and Analytical solution
figure
plot(ite3,wc3*1e+6,'LineWidth',1.2,'Color','r')
hold on
plot(ite1,wc1*1e+6,'LineWidth',1.2,'Color','b')
xlabel('Number of iterations')
ylabel('Membrane deflection at r = 0 (\mum)')
title('Non-uniform force Convergence Study')
legend('Averaged Analytical','Non-averaged FEM')
grid minor
set(gca,'FontSize',12)
set(gcf,'position',[100 150 800 400]);

%% Compare FEM 1D and COMSOL
mP = "N";
[w1,wuni1,ite1,wc1,rmp1,p1] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);
[w3,wuni3,ite3,wc3,rmp3] = analytical_1D(rm,nn,rb,E,d,T,mP);

% Load the COMSOL Model
model = mphload('2D_axisymm_st_rm_rb_diff.mph');

Rmem = int2str(rm*1e+3)+ "[mm]";
G = int2str((rm-rb)*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('G', G, 'Slit Gap between backplate and membrane');
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
w2unidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6]);
w2uniEf = mpheval(model,"es.Ez",'edim','boundary','selection',[3 6]);

% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);
model.component('comp1').common('sym1').active(true);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2nunidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6]);
w2nuniEf = mpheval(model,"es.Ez",'edim','boundary','selection',[3 6]);
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
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
legend('FEM Uniform','COMSOL Uniform')
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

subplot(2,2,2)
plot(rmp1*1e+3,wuni1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(rmp3*1e+3,wuni3*1e+6,'--','LineWidth',1.2,'Color','r')
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
legend('FEM Uniform','Analytical Uniform')
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Non-Uniform force
subplot(2,2,3)
plot(rmp1*1e+3,w1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(w2nunidisp.p(1,:)*1e+3,w2nunidisp.d1*1e+6,'--','LineWidth',1.2,'Color','r')
hold on
plot(rmp3*1e+3,wuni3*1e+6,'LineWidth',1.2,'Color','g')
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
legend('FEM Nonuniform','COMSOL Nonuniform','Analytical Uniform')
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Uniform Electric field
subplot(2,2,4)
plot(w2uniEf.p(1,:)*1e+3,w2uniEf.d1,'LineWidth',1.2,'Color','b')
hold on
plot(w2nuniEf.p(1,:)*1e+3,w2nuniEf.d1,'--','LineWidth',1.2,'Color','r')
xlabel('Radius (mm)')
ylabel('Electric Field (V/m)')
title('COMSOL')
legend('Uniform','Nonuniform')
grid minor
set(gca,'FontSize',12)
%set(gca,'position',[0.35 0.11 0.33 0.34])


% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Vs COMSOL 2D','t',supAxes);
set(gcf,'position',[50 50 870 680]);