clc
clear
close all

rm = 10e-3;      % radius of the membrane
rb = rm - 2e-3;         % radius of the backplate
E = 200;         % DC polarization voltage
d = 23.5*1e-6;   % distance between membrane and backplate
T = 3000;        % Tension of the membrane
bhr = 0.0005;    % Backplate hole radius
bhtheta = 30;    % Angle of separation between backplate holes
Nbh = 12;
thetaS = 15;     % Relative shift between two rings
bhrings = 4e-3;  % Distance from the hole to center of back plate
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
func = "VI";

%% Convergence study

cnt = 1;
for nn = 10:5:40
    
    tic;
    
    % Update the mesh each iteration
    hmax = rm/nn;    % Maximum length of the element
    mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);
    [p,e,t] = meshToPet(mesh);
    
    % Get the A and B matrices
    [w,~,~,xcirc,ycirc,~,~] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                    bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func);    
    t = toc;
    
    % Log the deflection at r = 0 and time
    wnuni(cnt) = mean(w);
    wmax(cnt) = max(abs(w));
    eta(cnt) = t;
    elelen(cnt) = "1/" + int2str(nn);
    cnt = cnt + 1;
end
%% Plot the convergence parameters

figure
subplot(1,2,1)
Nbh = 360/bhtheta;
for ii = 1:length(bhrings)
    for jj = 1:Nbh
        plot(xcirc(jj,:,ii)*1e+3,ycirc(jj,:,ii)*1e+3,'LineWidth',1.2,'Color','r')
        hold on
    end
end
xlabel('Backplate x direction (mm)')
ylabel('Backplate y direction (mm)')
xlim([-rb rb]*1e+3)
ylim([-rb rb]*1e+3)
title('Backplate Holes')
set(gca,'FontSize',12)

elevec = ["r_m/"+10 "r_m/"+20 "r_m/"+30 "r_m/"+40];
subplot(1,2,2)
yyaxis left
plot(wnuni*1e+6,'LineWidth',1.2)
ylabel('Mean Deflection (\mum)')
yyaxis right
plot(eta,'--','LineWidth',1.2)
ylabel('ETA (seconds)')
xlabel('Length of Element (Triangular Linear)')
title('Mesh Convergence Study')
xticks(1:2:cnt)
xticklabels(elevec)
grid on
set(gca,'FontSize',12)
set(gcf,'position',[100 150 900 400]);