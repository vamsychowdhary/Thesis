clc
clear
close all

rm = 10e-3;      % radius of the membrane
rb = rm - 2e-3;         % radius of the backplate
E = 200;         % DC polarization voltage
d = 23.5*1e-6;   % distance between membrane and backplate
T = 3000;        % Tension of the membrane
hmax = rm/30;    % Maximum length of the element
bhr = [0.0005 0.0005];    % Backplate hole radius
bhtheta = [30 30];    % Angle of separation between backplate holes
thetaS = 15;     % Relative shift between two rings
Nbh = [12 12];
bhrings = [4.5 6]*1e-3; % Distance from the hole to center of back plate
                 % A value more than one indicates concentric rings of
                 % holes.
bh = 1;          % Indicates if the backplate has holes or not. 0 means no holes and
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

func = "VI";
[w,~,~,xcirc,ycirc,~,~] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                   bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func);

%% Plot the membrane deflection
color = ["r" "b"];
figure
subplot(1,2,1)
for ii = 1:length(bhrings)
    for jj = 1:Nbh(ii)
        plot(xcirc(jj,:,ii)*1e+3,ycirc(jj,:,ii)*1e+3,'Color',color(ii),'LineWidth',1.2)
        hold on
    end
end
xlabel('Backplate x direction (mm)')
ylabel('Backplate y direction (mm)')
xlim([-rb rb]*1e+3)
ylim([-rb rb]*1e+3)
title('Backplate Holes')
set(gca,'FontSize',12)

subplot(1,2,2)
scatter3(p(1,:)*1e+3,p(2,:)*1e+3,w*1e+6,1)
xlabel('x direction (mm)')
ylabel('y direction (mm)')
zlabel('Deflection (\mum)')
title('FEM 2D')
zlim([-(0.5+max(abs(w*1e+6))) 0])
view(90,0)
set(gca,'FontSize',12)
set(gcf,'position',[100 150 700 350]);