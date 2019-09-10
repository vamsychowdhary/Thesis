function [w,wuni,ite,wc,rmp] = analytical_1D(rm,nn,rb,E,d,T,mP)

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
%
% Output Parameters
% w     -  Membrane deflection due to nonuniform pressure
% wuni -  Membrane deflection due to uniform pressure
% ite    -  Number of iterations (vector)
% wc    -  Membrane deflection at center for each ite
% rmp  -  Nodal coordinates

e0 = 8.85*1e-12; % permittivity of air

% Get the node points
rmp = linspace(0,rm,nn)'; % number of points on the membrane radial axis

% Calculate deflection
if (rb == rm)
        
    % Calculate the pressure due to electrostatic force
    xdorg = d*ones(nn,1);
    p = E^2*e0./(2*xdorg.^2);
    
    % Calculate the deflection
    w0 = p*rm^2/(4*T);
    wuni = -w0.*(1-(rmp/rm).^2);
    wtemp = wuni;
    
    % Now iterate using the non-uniform electrostatic force until the
    % deflection does not change
    cont = 't';
    numite = 0;
    while(cont == 't')
        
        % Calculate the new distance between membrane and backplate
        xd = xdorg+wtemp;
        
        % Calculate the average pressure as w0 is independant of the
        % position on the membrane
        p = E^2*e0./(2*xd.^2);
        if (mP == "Y")
            avgP = mean(p);
            p = (p./p)*avgP;
        end

        % Calculate the deflection
        w0 = p*rm^2/(4*T);
        w = -w0.*(1-(rmp/rm).^2);

        % compare the current deflection with previous deflection
        diffP = w(1)/wtemp(1);
        if (diffP > 1 && abs(w(1)) < d)
             wtemp = w;
             numite = numite + 1;
             ite(numite) = numite;
             wc(numite) = w(1);
        else
            w = wtemp;
            cont = 'f';
        end

    end

else
    
    % Find the points on backplate and slit gap
    ind = find(rmp <= rb);
    rbp = rmp(ind);
    nbp = length(rbp);
    clear ind
    ind = find(rmp > rb);
    rsp = rmp(ind);
    
    % Calculate the pressure due to electrostatic force
    xdorg = d*ones(nbp,1);
    p = E^2*e0./(2*xdorg.^2);
    avgP = mean(p);
    
    % Deflection for 0<r<rb
    wb = (avgP/(4*T))*(rbp.^2 + rb^2*(2*log(rb/rm) - 1));
    
    % Deflection for rb<r<rm
    ws = (avgP*rb^2/(2*T))*log(rsp/rm);
    
    % Total deflection due to uniform Electrostatic pressure
    wuni = [wb' ws'];
    wbtemp = wb;
    
    % Now iterate using the non-uniform electrostatic force until the
    % deflection does not change
    cont = 't';
    numite = 0;
    while(cont == 't')
        
        % Calculate the new distance between membrane and backplate
        xd = xdorg+wbtemp;
        
        % Calculate the average pressure as w0 is independant of the
        % position on the membrane
        p = E^2*e0./(2*xd.^2);
        avgP = mean(p);
        
        % Deflection for 0<r<rb
        wb = (avgP/(4*T))*(rbp.^2 + rb^2*(2*log(rb/rm) - 1));

        % Deflection for rb<r<rm
        ws = (avgP*rb^2/(2*T))*log(rsp/rm);
        
        % Total deflection due to uniform Electrostatic pressure
        w = [wb' ws'];
        
        % compare the current deflection with previous deflection
        diffP = w(1)/wbtemp(1);
        if (diffP > 1)
             wbtemp = wb;
             wtemp = w;
             numite = numite + 1;
             ite(numite) = numite;
             wc(numite) = w(1);
        else
            w = wtemp;
            cont = 'f';
        end

    end
    
end

end