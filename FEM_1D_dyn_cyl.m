function [w,rmp,p,e,edB] = FEM_1D_dyn_cyl(rm,nn,rb,E,d,T,func,f,pin,surfd,wini)

e0 = 8.85*1e-12; % permittivity of air
cm = T/surfd;
w = 2*pi*f;
K = w^2/cm;

% Get the node points
rmp = linspace(0,rm,nn); % number of points on the membrane radial axis

% Get the A and B matrices
if (func == "VA")
    [A,B] = FEM_1D_Linear_A_B_Global(rmp,nn);
elseif (func == "VI")
    [A,B] = FEMmemLIN(rmp,nn);
end

% input pressure
p = pin*ones(nn,1);

% Boundary condition - Set the deflection on the edges to zero
A(nn,nn-1:nn) = [0 1];
B(nn,nn-1:nn) = [0 0];

w = zeros(2,nn);

% Calculate the deflection without static
w(2,:) = (((A+K*B)*T)\(B*p))';

% Calculate the deflection with static
w(1,:) =  w(2,:)+ wini';

% Calculate the sensitivity with static deflection
ind = find(rmp <= rb);
Cs = e0*2*pi*rmp(2)*sum(rmp(ind)'./(d+wini(ind)));
Q = Cs*E;
Cd = e0*2*pi*rmp(2)*sum(rmp(ind)'./(d+w(1,ind)'));
e(1) = E - Q/Cd;
edB(1) = 20*log10(abs(e(1)));

% Calculate the sensitivity without static deflection
dn = d*ones(max(ind),1);
Cs = e0*2*pi*rmp(2)*sum(rmp(ind)'./dn);
Q = Cs*E;
Cd = e0*2*pi*rmp(2)*sum(rmp(ind)'./(d+w(2,ind)'));
e(2) = E - Q/Cd;
edB(2) = 20*log10(abs(e(2)));
end