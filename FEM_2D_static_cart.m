function [w,wuni,ite,xcirc,ycirc,A,B] = FEM_2D_static_cart(rm,rb,E,d,T,...
                                        bh,bhr,bhtheta,thetaS,bhrings,th,Nbh,p,e,t,func)
e0 = 8.85*1e-12; % permittivity of air

% Get the A and B matrices
if (func == "VA")
    [A,B] = FEM_2D_TRI_Linear_A_B(p,t);
elseif (func == "VI")
    [A,B] = FEMmem2D(p,t);
end

% Boundary condition - Set the deflection on the edges to zero
nn = size(p,2);
rd = sqrt(p(1,:).^2+p(2,:).^2);
ind = e(1,:); % Geometric edges from Mesh
for ii = 1:length(ind)
    A(ind(ii),:) = [zeros(1,ind(ii)-1) 1 zeros(1,nn-ind(ii))];
    B(ind(ii),:) = zeros(1,nn);
end

% Calculate the pressure due to electrostatic force
xdorg = d*ones(1,nn);
pe = E^2*e0./(2*xdorg.^2);

% check membrane and back plate radius
clear ind
if (rb < rm)
    ind = find(rb < abs(rd));
    pe(ind) = 0;
end

% If the back plate has holes make pe zero at those points
xcirc = zeros(max(Nbh),length(th),length(bhrings));
ycirc = zeros(max(Nbh),length(th),length(bhrings));
indbh = 0;

% check if holes present on the backplate
if(bh == 1)
    % Loop through each annular ring of holes
     theta = 0;
    for ii = 1:length(bhrings)
        % Loop through all the holes in an annular ring
        for jj = 1:Nbh(ii)
            xbh = bhrings(ii)*cosd(theta);
            ybh = bhrings(ii)*sind(theta);
            bhrd = sqrt((p(1,:)-xbh).^2+(p(2,:)-ybh).^2);
            indn = find(bhrd < bhr(ii));
            pe(indn) = 0;
            indbh = [indbh indn];
            theta = theta + bhtheta(ii);
            
            % Break if the hole is at center
            if (bhrings(ii) == 0)
                
                % Only for plotting
                xcirc(jj,:,ii) = bhr(ii)*cos(th) ;
                ycirc(jj,:,ii) = bhr(ii)*sin(th) ;
                xcirc(jj+1:end,:,ii) = 0;
                ycirc(jj+1:end,:,ii) = 0;                
                break;
            end
            
            % Only for plotting
            xcirc(jj,:,ii) = bhr(ii)*cos(th) + xbh;
            ycirc(jj,:,ii) = bhr(ii)*sin(th) + ybh;
        end
        theta = theta + thetaS;
    end
end

% Calculate the deflection
wuni = (A*T)\(B*pe');
wtemp = wuni;

% Now iterate using the non-uniform electrostatic force until the
% deflection does not change
cont = 't';
numite = 0;
indC = find(wuni == min(wuni));
while(cont == 't')
    
    xd = xdorg+wtemp';
    pe = E^2*e0./(2*xd.^2);
    
    % Check membrane and back plate radius
    if (rb < rm )
        pe(ind) = 0;
    end
    
    % Check for holes on backplate
    if (bh == 1)
        pe(indbh(2:end)) = 0;
    end
    
    % Calculate the deflection
    w = (A*T)\(B*pe');
    
    % compare the current deflection with previous deflection
    %diffP = abs(mean(w))/abs(mean(wtemp));
    diffP = w(indC)/wtemp(indC);
    if (diffP > 1 && w(indC) < d)
        wtemp = w;
         numite = numite + 1;
         ite(numite) = numite;
         wc(numite) = w(indC);
    else
        if (w(indC) > d)
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