function [w,wuni,ite,wc,rmp,p] = FEM_1D_static_cyl(rm,nn,rb,E,d,T,mP,rw,rp,func)

% This function calculates the static deflection of a membrane for a 1D
% axisymmetrical geometry.
%
% Input Parameters
% rm  - Radius of the membrane
% nn  - Number of nodes
% rb   - Radius of the backplate
% E    - Polarization voltage
% d    - Physical distance between membrane and backplate
% T    - Tension of the membrane
% mP  - Apply mean pressure - Y or N
% rw   - Width of the annular ring
% rp    - Position of the annular ring from the center
% func - Function to use. VA - Vamsy and VI - Vicente
%
% Output Parameters
% w     -  Membrane deflection due to nonuniform pressure
% wuni -  Membrane deflection due to uniform pressure
% ite    -  Number of iterations (vector)
% wc    -  Membrane deflection at center for each ite
% rmp  -  Nodal coordinates
% p      -  Electrostatic pressure

e0 = 8.85*1e-12; % permittivity of air

% Get the node points
rmp = linspace(0,rm,nn); % number of points on the membrane radial axis

% Get the A and B matrices
if (func == "VA")
    [A,B] = FEM_1D_Linear_A_B_Global(rmp,nn);
elseif (func == "VI")
    [A,B] = FEMmemLIN(rmp,nn);
end

% Calculate the pressure due to electrostatic force
xdorg = d*ones(nn,1);
p = E^2*e0./(2*xdorg.^2);

% check membrane and back plate radius
if (rb < rm)
    ind = find(rmp > rb);
    nbp = nn-length(ind);
    p(ind) = 0;
end

% Find the points in the ring width and set the electrostatic pressure to
% zero
if (rw ~= 0 && rp ~= 0)
    indR = find((rmp > rp) & (rmp < (rp+rw)));
    p(indR) = 0;
else
    indR = 0;
end

% Boundary condition - Set the deflection on the edges to zero
A(nn,nn-1:nn) = [0 1];
B(nn,nn-1:nn) = [0 0];

% Calculate the deflection
wuni = (A*T)\(B*p);
wtemp = wuni;

% Now iterate using the non-uniform electrostatic force until the
% deflection does not change
cont = 't';
numite = 0;
while(cont == 't')
    
    % Update the distance and calculate the new electrostatic pressure
    xd = xdorg+wtemp;
    p = E^2*e0./(2*xd.^2);
    
    % check membrane and back plate radius
    if (rb == rm && mP == "Y")
        avgP = mean(p);
        p = (p./p)*avgP;
    elseif (rb < rm && mP == "Y")
        p(ind) = 0;
        avgP = sum(p(1:nbp))/nbp;
        p(1:nbp) = (p(1:nbp)./p(1:nbp))*avgP;
    elseif (rb < rm && mP == "N")
        p(ind) = 0;
    end
    
    % Check for the annular ring
    if (indR ~= 0)
        p(indR) = 0;
    end
    
    % Calculate the deflection
    w = (A*T)\(B*p);
    
    % compare the current deflection with previous deflection
    diffP = w(1)/wtemp(1);
    if (diffP > 1 && abs(w(1)) < d)
         wtemp = w;
         numite = numite + 1;
         ite(numite) = numite;
         wc(numite) = w(1);
    else
        % If the deflection is more than the physical separaration -
        % membrane colapsed to backplate
        if (abs(w(1)) > d)
           msg = "Membrane collapsed to backplate. Please consider"+...
                " increasing the static Tension or decreasing the Polarization voltage"+...
                " or increasing the air gap width or decreasing the membrane radius";
           disp(msg)
        end
        w = wtemp;
        cont = 'f';            
    end
       
end

end