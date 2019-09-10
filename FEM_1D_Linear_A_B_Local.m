function [A,B] = FEM_1D_Linear_A_B_Local(rmp,nn)

% This function calculates the A and B matrices for given Linear element nodes
% using Local coordinates
% 
% Input Parameters
% rmp - Nodal coordinates of the linear elements
% nn   - Number of nodes
% 
% Output Parameters
% A     - Assembled A matrix
% B     - Assembled B matrix

% Initialize the matrices
A = zeros(nn,nn);
B = zeros(nn,nn);
Imat = [1 -1;-1 1];
B1 = [8/3 4/3;4/3 8/3];
B2 = [-4/3 0;0 4/3];

% Loop through each element
for ii = 1:nn-1
    
    % Constants for a given element
    h1 = rmp(ii+1)-rmp(ii);
    h2 = rmp(ii+1)+rmp(ii);
    
    % Compute A matrix and assemble it
    A(ii:ii+1,ii:ii+1) = A(ii:ii+1,ii:ii+1)-(h2/(2*h1))*Imat;
    
    % Compute B matrix and assemble it
    Bele = (B1*h1*h2 + B2*h1^2)/16;
    B(ii:ii+1,ii:ii+1) = B(ii:ii+1,ii:ii+1) + Bele;
end

end