clc
clear
close all

rm = 10e-3;      % radius of the membrane
rb = rm - 2e-3;         % radius of the backplate
E = 200;         % DC polarization voltage
d = 23.5*1e-6;   % distance between membrane and backplate
T = 3000;        % Tension of the membrane
hmax = rm/30;    % Maximum length of the element
bhr = 0.0005;    % Backplate hole radius
bhtheta = 30;    % Angle of separation between backplate holes
bhrings = 4e-3; % Distance from the hole to center of back plate
                 % A value more than one indicates concentric rings of
                 % holes.
bh = 1;          % Indicates if the backplate has holes or not. 0 means no holes and
                 % 1 means yes.
th = 0:pi/50:2*pi; % Number of points on the cricle for plotting
Nbh = 12;
thetaS = 0;

% Using MATLAB PDE Toolbox to generate mesh
gd = [4;0;0;rm;rm;0];
ns = [69;49];
sf = 'E1';
g = decsg(gd,sf,ns);
model = createpde;
geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);
[p,e,t] = meshToPet(mesh);

func = "VI";
[w1,wuni1,ite1,~,~,A1,B1] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                      bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func);

% Load the COMSOL Model
model = mphload('3D_model_st_rm_rb_diff_holes.mph');

Rmem = int2str(rm*1e+3)+ "[mm]";
Vpol = int2str(E)+"[V]";
Tm0 = int2str(T)+"[N/m]";
Hm = num2str(d*1e+6,"%3.1f")+"[um]";
Rhole = num2str(bhr*1e+3,"%2.1f")+"[mm]";
Rang = num2str(bhtheta*pi/180,"%5.4f");
Rdist = int2str(bhrings*1e+3)+"[mm]";

% Update model parameters accordingly
model.param.set('Rmem', Rmem, 'Radius of membrane');
model.param.set('Vpol',Vpol,'Polarization voltage');
model.param.set('Tm0',Tm0,'Membrane static tension');
model.param.set('Hm',Hm,'Air gap width');
model.param.set('Rhole',Rhole,'Radius of hole on the backplate');
model.param.set('Hm',Hm,'Angle of separation between two holes');
model.param.set('Rdist',Rdist,'Distance of hole from center');


% Enable the moving mesh components to get deflection due to non-uniform force
model.component('comp1').common('free1').active(true);
model.component('comp1').common('fix1').active(true);
model.component('comp1').common('disp1').active(true);

% Run the model
model.study("std1").run

% Load the displacement z component and e field to a variable
w2nunidisp = mpheval(model,"w",'dataset','dset2','edim','boundary','selection',4);
%% Plot the membrane deflection
figure
subplot(1,2,1)
scatter3(p(1,:)*1e+3,p(2,:)*1e+3,w1*1e+6,2)
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
title('FEM 2D')
zlim([-(0.5+max(abs(w1*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)
subplot(1,2,2)
scatter3(w2nunidisp.p(1,:)*1e+3,w2nunidisp.p(2,:)*1e+3,w2nunidisp.d1*1e+6,2)
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
title('COMSOL 3D')
zlim([-(0.5+max(abs(w1*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)
set(gcf,'position',[100 150 700 350]);
%% Plot the back plate holes
Nbh = 360/bhtheta;
figure
for ii = 1:length(bhrings)
    for jj = 1:Nbh
        plot(xcirc2(jj,:,ii),ycirc2(jj,:,ii),'Color','b')
        hold on
    end
end