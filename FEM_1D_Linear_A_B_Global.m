function [A,B] = FEM_1D_Linear_A_B_Global(rmp,nn)

% This function calculates the A and B matrices for given Linear element nodes
% using Global coordinates
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
Bele = zeros(2,2);
Imat = [1 -1;-1 1];

% Loop through each element
for ii = 1:nn-1
    
    % Load the nodes to temp vars
    a = rmp(ii);
    b = rmp(ii+1);
    
    % Compute A matrix and assemble it
    A(ii:ii+1,ii:ii+1) = A(ii:ii+1,ii:ii+1)-(a+b)*Imat;
    
    % Compute each value of Bmat for a given element
    Bele(1,1) = (b^4 - 6*(a*b)^2 + 8*a^3*b - 3*a^4)/6;
    Bele(1,2) = (b^4 - 2*a*b^3 + 2*a^3*b - a^4)/6;
    Bele(2,1) = Bele(1,2);
    Bele(2,2) = (3*b^4 - 8*a*b^3 + 6*(a*b)^2 - a^4)/6;
    
    % Assemble B matrix
    B(ii:ii+1,ii:ii+1) = B(ii:ii+1,ii:ii+1) + Bele/(b-a);
end

end