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
rw = 2e-3;         % Width of the concentric ring
rp = 2e-3;         % Position of the concentric ring from the center
func = "VA";       % Use Vamsy's function

%% Compare FEM 1D and COMSOL
mP = "N";
[w1,wuni1,ite1,wc1,rmp1,p1] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);

% Load the COMSOL Model
model = mphload('2D_axisymm_st_rm_rb_diff_ring.mph');

Rmem = int2str(rm*1e+3)+ "[mm]";
G = int2str((rm-rb)*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";
Rw = int2str(rw*1e+3)+ "[mm]";
Rp = int2str(rp*1e+3)+ "[mm]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('G', G, 'Slit Gap between backplate and membrane');
model.param.set('Vpol',Vpol,'Polarization voltage');
model.param.set('Tm0',Tm0,'Membrane static tension');
model.param.set('Hm',Hm,'Air gap width');
model.param.set('rw',Rw,'Ring width');
model.param.set('rp',Rp,'Ring position from the center');

% Disable the moving mesh components to get deflection due to uniform force
model.component('comp1').common('free1').active(false);
model.component('comp1').common('fix1').active(false);
model.component('comp1').common('disp1').active(false);
model.component('comp1').common('sym1').active(false);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2unidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6 9 12]);
w2uniEf = mpheval(model,"es.Ez",'edim','boundary','selection',[3 6 9 12]);

% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);
model.component('comp1').common('sym1').active(true);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2nunidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6 9 12]);
w2nuniEf = mpheval(model,"es.Ez",'edim','boundary','selection',[3 6 9 12]);
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
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend('FEM Nonuniform','COMSOL Nonuniform');
lgd.FontSize = 10;
ylim([maxy 0])
grid minor
set(gca,'FontSize',12)

% Compare Uniform and Non-uniform Electric field
subplot(2,2,3)
plot(w2uniEf.p(1,:)*1e+3,w2uniEf.d1,'LineWidth',1.2,'Color','b')
hold on
plot(w2nuniEf.p(1,:)*1e+3,w2nuniEf.d1,'--','LineWidth',1.2,'Color','r')
xlabel('Radius (mm)')
ylabel('Electric Field (V/m)')
title('COMSOL')
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

%% FEM 1D - Compare the deflections for same ring width but different distances from center of membrane

mP = "N";
rwn = 1e-3;         % Width of the concentric ring
rp = [1 1.5 3 4 4.5 6]*1e-3;         % Position of the concentric ring from the center
cnt = 0;
for ii = 1:length(rp)
    cnt = cnt+1;
    [w(cnt,:),~,~,~,rmp(cnt,:),~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rwn,rp(ii),func);
end

%% FEM 1D - Compare the deflections for varying ring widths but same distance from center of membrane

mP = "N";
rw = [1 1.5 3 1 1.5 3]*1e-3;         % Width of the concentric ring. 
rpn = 1e-3;         % Position of the concentric ring from the center
rpn2 = 4e-3;
for ii = 1:length(rw)
    cnt = cnt+1;
    if (ii > 3)
        rpn = rpn2;
    end
    [w(cnt,:),~,~,~,rmp(cnt,:),~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw(ii),rpn,func);
end

%% Plot the results from the above two sections

figure
subplot(1,2,1)
for ii = 1:length(rp)
    plot(rmp(ii,:)*1e+3,w(ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rp = "+int2str(rp(ii)*1e+3)+"mm, Rw = "+int2str(rwn*1e+3)+"mm";
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
for ii = 1:length(rw)
    plot(rmp(length(rp)+ii,:)*1e+3,w(length(rp)+ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rw = "+int2str(rw(ii)*1e+3)+"mm Rp = "+int2str(rpn*1e+3)+"mm";
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.09 0.09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Nonuniform Pressure','t',supAxes);
set(gcf,'position',[50 50 800 350]);

%% Plot the case of changing ring position only
clear lgdv
e0 = 8.85*1e-12; % permittivity of air
figure
subplot(1,2,1)
for ii = 1:3
    plot(rmp(ii,:)*1e+3,w(ii,:)*1e+6,'LineWidth',1.2)

    % Calculate the mean electrostatic pressure
    p = E^2*e0./(2*(d+w(ii,:)).^2);
    indR = find((rmp(ii,:) > rp(ii)) & (rmp(ii,:) < (rp(ii)+rw)));
    p(indR) = 0;
    MP = mean(p);
%     lgdv(ii) = "Rp = "+num2str(rp(ii)*1e+3,'%2.1f')+"mm, Rw = "+num2str(rwn*1e+3,'%2.1f')+"mm, MEP: " ...
%                         + num2str(MP,'%4.2f') ;
    lgdv(ii) = "Rp = "+num2str(rp(ii)*1e+3,'%2.1f')+"mm, Rw = "+num2str(rwn*1e+3,'%2.1f')+"mm" ;
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
for ii = 4:6
    plot(rmp(ii,:)*1e+3,w(ii,:)*1e+6,'LineWidth',1.2)
    
    % Calculate the mean electrostatic pressure
    p = E^2*e0./(2*(d+w(ii,:)).^2);
    indR = find((rmp(ii,:) > rp(ii)) & (rmp(ii,:) < (rp(ii)+rw)));
    p(indR) = 0;
    MP = mean(p);
%     lgdv(ii) = "Rp = "+num2str(rp(ii)*1e+3,'%2.1f')+"mm, Rw = "+num2str(rwn*1e+3,'%2.1f')+"mm, MEP: " ...
%                         + num2str(MP,'%4.2f') ;
    lgdv(ii) = "Rp = "+num2str(rp(ii)*1e+3,'%2.1f')+"mm, Rw = "+num2str(rwn*1e+3,'%2.1f')+"mm";
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv(4:6));
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.09 0.09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Nonuniform Pressure','t',supAxes);
set(gcf,'position',[50 50 800 350]);

%% Plot the case of changing ring width only
rpn = 1e-3;
clear lgdv
figure
subplot(1,2,1)
for ii = 1:3
    plot(rmp(ii,:)*1e+3,w(ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rw = "+num2str(rw(ii)*1e+3,'%2.1f')+"mm, Rp = "+num2str(rpn*1e+3,'%2.1f')+"mm";
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
for ii = 4:6
    plot(rmp(ii,:)*1e+3,w(ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rw = "+num2str(rw(ii)*1e+3,'%2.1f')+"mm, Rp = "+num2str(rpn2*1e+3,'%2.1f')+"mm";
    hold on
end
hold off
xlabel('Radius (mm)')
ylabel('Deflection (\mum)')
lgd = legend(lgdv(4:6));
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.09 0.09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Nonuniform Pressure','t',supAxes);
set(gcf,'position',[50 50 800 350]);

%% COMSOL 2D - Compare the deflections for same ring width but different distances from center of membrane

rwn = 1e-3;         % Width of the concentric ring
rp = [2 4 6]*1e-3;         % Position of the concentric ring from the center
cnt = 0;

% Load the COMSOL Model
model = mphload('2D_axisymm_st_rm_rb_diff_ring.mph');

% Set the model parameters
Rmem = int2str(rm*1e+3)+ "[mm]";
G = int2str((rm-rb)*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";
Rw = int2str(rwn*1e+3)+ "[mm]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('G', G, 'Slit Gap between backplate and membrane');
model.param.set('Vpol',Vpol,'Polarization voltage');
model.param.set('Tm0',Tm0,'Membrane static tension');
model.param.set('Hm',Hm,'Air gap width');
model.param.set('rw',Rw,'Ring width');

% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);
model.component('comp1').common('sym1').active(true);

for ii = 1:length(rp)
    cnt = cnt+1;
    
    % Update the ring position in each iteration
    Rp = int2str(rp(ii)*1e+3)+ "[mm]";
    model.param.set('rp',Rp,'Ring position from the center');
    
    % Run the model
    model.study("std1").run
    
    % Load the displacement z component and e field to a variable
    w = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6 9 12]);
    wC(cnt,:) = w.d1;
    rmpC(cnt,:) = w.p(1,:);
end

%% COMSOL 2D - Compare the deflections for different ring widths but same distance from center of membrane

rw = [1 2 4]*1e-3;         % Width of the concentric ring
rpn = 2e-3;         % Position of the concentric ring from the center

% Set the model parameters
Rmem = int2str(rm*1e+3)+ "[mm]";
G = int2str((rm-rb)*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";
Rp = int2str(rpn*1e+3)+ "[mm]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('G', G, 'Slit Gap between backplate and membrane');
model.param.set('Vpol',Vpol,'Polarization voltage');
model.param.set('Tm0',Tm0,'Membrane static tension');
model.param.set('Hm',Hm,'Air gap width');
model.param.set('rw',Rp,'Ring position from center');

for ii = 1:length(rw)
    cnt = cnt+1;
    
    % Update the ring position in each iteration
    Rw = int2str(rw(ii)*1e+3)+ "[mm]";
    model.param.set('rw',Rw,'Ring width');
    
    % Run the model
    model.study("std1").run
    
    % Load the displacement z component and e field to a variable
    w = mpheval(model,"spatial.dz",'edim','boundary','selection',[3 6 9 12]);
    wC(cnt,:) = w.d1;
    rmpC(cnt,:) = w.p(1,:);
end

%% Plot the results from the above two sections

figure
subplot(1,2,1)
for ii = 1:length(rp)
    plot(rmpC(ii,:)*1e+3,wC(ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rp = "+int2str(rp(ii)*1e+3)+"mm, Rw = "+int2str(rwn*1e+3)+"mm";
    hold on
end
hold off
ylabel('Membrane deflection (\mum)')
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

subplot(1,2,2)
for ii = 1:length(rw)
    plot(rmpC(length(rp)+ii,:)*1e+3,wC(length(rp)+ii,:)*1e+6,'LineWidth',1.2)
    lgdv(ii) = "Rw = "+int2str(rw(ii)*1e+3)+"mm Rp = "+int2str(rpn*1e+3)+"mm";
    hold on
end
hold off
xlabel('Membrane radial direction (mm)','position',[-1.5 -2.67 -1])
lgd = legend(lgdv);
lgd.FontSize = 10;
grid minor
set(gca,'FontSize',12)

% Common X and Y labels
% 1st value in supAxes specifies the distance of Ylabel from the plot
% 2nd value in supAxes specifies the distance of Xlabel from the plot
supAxes = [.09 0.09 0.84 0.84];
% suplabel('Membrane radial direction (mm)','x',supAxes);
% suplabel('Membrane deflection (\mum)','y',supAxes);
suplabel('FEM 1D Non-uniform force','t',supAxes);
set(gcf,'position',[50 50 800 350]);