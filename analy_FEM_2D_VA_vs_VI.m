clc
clear
close all

rm = 10e-3;      % radius of the membrane
rb = rm;         % radius of the backplate
E = 200;         % DC polarization voltage
d = 23.5*1e-6;   % distance between membrane and backplate
T = 3000;        % Tension of the membrane
hmax = rm/30;    % Maximum length of the element
bhr = 0.0005;    % Backplate hole radius
bhtheta = 30;    % Angle of separation between backplate holes
thetaS = 15;     % Relative shift between two rings
Nbh = [1 6 12];    % No holes in each ring
bhrings = [0 0.002 0.003 0.004 0.005 0.006]; % Distance from the hole to center of back plate
                 % A value more than one indicates concentric rings of
                 % holes.
bh = 0;          % Indicates if the backplate has holes or not. 0 means no holes and
                 % 1 means yes.
th = 0:pi/50:2*pi; % Number of points on the cricle for plotting

% Using MATLAB PDE Toolbox to generate mesh
gd = [4;0;0;rm;rm;0];
ns = [69;49];
sf = 'E1';
g = decsg(gd,sf,ns);
model = createpde;
geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);
[p,e,t] = meshToPet(mesh);

func = "VA";
[w1,wuni1,ite1,xcirc,ycirc,A1,B1] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                      bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func);


func = "VI";
[w2,wuni2,ite2,~,~,A2,B2] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                      bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func);

hmax = rm/90;    % Maximum length of the element
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);
[p1,e1,t1] = meshToPet(mesh);

func = "VA";
[w3,wuni3,ite3,~,~,A3,B3] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                      bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p1,e1,t1,func);
%%  Find the deflection on the cross section
cnt = 0;
for ii = 1:size(p,2)
    if (p(1,ii) >= 0 && p(2,ii) == 0)
        cnt = cnt+1;
        ind(cnt) = ii; 
    end
end
plot(p(1,ind),w1(ind))
%% COMSOL

% Load the COMSOL Model
model = mphload('3D_model_st_rm_rb_same.mph');

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

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2unidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',4);
w2uniEf = mpheval(model,"es.Ez",'edim','boundary','selection',4);

% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2nunidisp = mpheval(model,"spatial.dz",'edim','boundary','selection',4);
w2nuniEf = mpheval(model,"es.Ez",'edim','boundary','selection',4);
%%  Call the FEM 1D
func = "VA";
nn = 2000;
mP = "N";   % Apply a mean pressure
rw = 0;
rp = 0;
[w4,~,~,~,rmp4,~] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func);

%% Plot the membrane deflection
figure
subplot(2,2,1)
scatter3(p(1,:)*1e+3,p(2,:)*1e+3,w1*1e+6,1)
hold on
scatter3(rmp4*1e+3,rmp4*1e+3,w4*1e+6,1,'MarkerEdgeColor','k')
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
legend('2D','1D')
title('FEM')
zlim([-(0.5+max(abs(w2*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)

subplot(2,2,2)
scatter3(p1(1,:)*1e+3,p1(2,:)*1e+3,w3*1e+6,1)
hold on
scatter3(rmp4*1e+3,rmp4*1e+3,w4*1e+6,1,'MarkerEdgeColor','k')
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
legend('2D','1D')
title('FEM')
zlim([-(0.5+max(abs(w2*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)

subplot(2,2,3)
scatter3(p(1,:)*1e+3,p(2,:)*1e+3,w2*1e+6,1,'MarkerEdgeColor','r')
hold on
scatter3(rmp4*1e+3,rmp4*1e+3,w4*1e+6,1,'MarkerEdgeColor','k')
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
title('FEM')
legend('2D','1D')
zlim([-(0.5+max(abs(w2*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)

% subplot(2,2,4)
% plot(rmp4*1e+3,w4*1e+6,'LineWidth',1.2,'Color','k')
% xlabel('Radius (mm)')
% ylabel('Deflection (\mum)')
% title('FEM 1D')
% set(gca,'FontSize',12)
set(gcf,'position',[100 100 870 680]);
supAxes = [.09 0.09 0.84 0.84];
%suplabel('FEM 2D','t',supAxes);
%% Plot the back plate holes
Nbh = 360/bhtheta;
figure
for ii = 1:length(bhrings)
    for jj = 1:Nbh
        plot(xcirc2(jj,:,ii),ycirc2(jj,:,ii),'Color','b')
        hold on
    end
end

%% Plot the mesh
pdemesh(model)