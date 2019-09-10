clc
clear
close all

rm = 10e-3;  % radius of the membrane
rb = rm - 2e-3;
nn = 3001;  % number of nodes
mP = "N";
rw = 0;
rp = 0;
func = "VA";

%% Different polarization voltages (E), gap widths (d) and Tension (T)

% Varying E
E1 = 200;            % DC voltage
T = 3000;           % Tension of the membrane
d = 23.5*1e-6;      % distance between membrane and backplate
[wE1,~,~,~,rmp1,~] = FEM_1D_static_cyl(rm,nn,rb,E1,d,T,mP,rw,rp,func);

E2 = 100;
[wE2,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E2,d,T,mP,rw,rp,func);

E3 = 50;
[wE3,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E3,d,T,mP,rw,rp,func);

% Varying d
E = 200;
d1 = 30*1e-6;
[wd1,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E,d1,T,mP,rw,rp,func);

d2 = 50*1e-6;
[wd2,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E,d2,T,mP,rw,rp,func);

% Varying T
d = 23.5*1e-6;
T1 = 2500;
[wT1,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T1,mP,rw,rp,func);

T2 = 2000;
[wT2,~,~,~,~,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T2,mP,rw,rp,func);

%% Plot the membrane deflection
figure
subplot(2,2,1)
plot(rmp1*1e+2,wE1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(rmp1*1e+2,wE2*1e+6,'--','LineWidth',1.2,'Color','r')
hold on
plot(rmp1*1e+2,wE3*1e+6,'LineWidth',1.2,'Color','g')
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
title('Varying Polarization Voltage (E)')
legend("Reference capsule",int2str(E2)+"V",int2str(E3)+"V")
grid minor
set(gca,'FontSize',12)

subplot(2,2,2)
plot(rmp1*1e+2,wE1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(rmp1*1e+2,wd1*1e+6,'--','LineWidth',1.2,'Color','r')
hold on
plot(rmp1*1e+2,wd2*1e+6,'LineWidth',1.2,'Color','g')
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
legend("Reference capsule",num2str(d1*1e+6,'%2.0f')+"\mum" ...
,num2str(d2*1e+6,'%2.0f')+"\mum")
title('Varying Air Gap Width (d)')
grid minor
set(gca,'FontSize',12)

subplot(2,2,3)
plot(rmp1*1e+2,wE1*1e+6,'LineWidth',1.2,'Color','b')
hold on
plot(rmp1*1e+2,wT1*1e+6,'--','LineWidth',1.2,'Color','r')
hold on
plot(rmp1*1e+2,wT2*1e+6,'LineWidth',1.2,'Color','g')
legend("Reference capsule",int2str(T1)+"N/m",int2str(T2)+"N/m")
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
title('Varying Tension (T)')
grid minor
set(gca,'FontSize',12)
set(gca,'position',[0.3325 0.13 0.334 0.341]);
set(gcf,'position',[100 100 870 680]);

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.12 .13 0.75 0.8];
%suplabel('Membrane radial direction (mm)','x',supAxes);
%suplabel('Membrane deflection (\mum)','y',supAxes);