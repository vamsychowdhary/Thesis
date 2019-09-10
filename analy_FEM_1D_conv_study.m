% In this script, the convergence study with respect to number of elements
% is studied for a 1D axisymmetric membrane by adjusting rb.

clc
clear
close all

e0 = 8.85*1e-12;      % Permittivity of air
rm = 10e-3;             % Radius of the membrane
E = [200 100];         % Polarization voltage
T = 3000;                % Static Tension of the membrane
d = 23.5*1e-6;         % Physical sepration between membrane and backplate
rb = rm - 2e-3;                  % Radius of the backplate
mP = "N";                % Apply mean pressure
rw = 0;                    % Width of the annular ring
rp = 0;                     % Position of the annular ring
func = "VA";              % function to use. VA - Vamsy, VI - Vicente

%% Convergence study

for jj = 1:length(E)
    ii = 1;
    for nn = 100:200:4000

        tic;

        % Get the A and B matrices
        [w,wuni,ite,wc,rmp,p] = FEM_1D_static_cyl(rm,nn,rb,E(jj),d,T,mP,rw,rp,func);

        t = toc;

        % Log the deflection at r = 0 and time
        wnuni(jj,ii) = w(1);
        eta(jj,ii) = t;
        nofele(ii) = nn-1;
        ii = ii + 1;
    end
end
%% Plot the convergence parameters

figure
subplot(1,2,1)
yyaxis left
plot(nofele,wnuni(1,:)*1e+6,'LineWidth',1)
ylabel('Membrane Deflection at rm = 0 (\mum)')
yyaxis right
plot(nofele,eta(1,:),'--','LineWidth',1)
ylabel('ETA (seconds)')
xlabel('Number of Elements (Linear)')
legend('E = 200V')
grid on
set(gca,'FontSize',12)

subplot(1,2,2)
yyaxis left
plot(nofele,wnuni(2,:)*1e+6,'LineWidth',1)
yyaxis right
plot(nofele,eta(2,:),'--','LineWidth',1)
xlabel('Number of Elements (Linear)')
title('Mesh Convergence Study')
legend('E = 100V')
grid on
set(gca,'FontSize',12)
set(gcf,'position',[100 150 900 370]);