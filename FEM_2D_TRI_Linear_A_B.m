function [A,B] = FEM_2D_TRI_Linear_A_B(p,t)

% Initialize the matrices
nn = size(p,2);
A = zeros(nn,nn);
B = zeros(nn,nn);
BtempC = [ 1/6   -1/12   -1/12 ;
          -1/12   1/3     1/4  ;
          -1/12   1/4     1/3] ;

% Loop through all the elements
for ii = 1:size(t,2)
    
    % locate the coordinates of each vertex
    l = t(1,ii);
    m = t(2,ii);
    n = t(3,ii);
    x1 = p(1,l);
    y1 = p(2,l);
    x2 = p(1,m);
    y2 = p(2,m);
    x3 = p(1,n);
    y3 = p(2,n);
    
    % Find the jacobian matrix and it's inverse
    J = [-x1+x2 -y1+y2;-x1+x3 -y1+y3];
    detJ = det(J);
    R = inv(J);
    
    % Compute the derivatives of the shape functions
    dwx = [-R(1,1)-R(1,2);R(1,1);R(1,2)];
    dwy = [-R(2,1)-R(2,2);R(2,1);R(2,2)];
    
    % Compute the conduction matrix and add to the corresponding global
    % elements
    Atemp = -(dwx*dwx' + dwy*dwy')*detJ;
    A(l,l) = A(l,l) + Atemp(1,1);
    A(l,m) = A(l,m) + Atemp(1,2);
    A(l,n) = A(l,n) + Atemp(1,3);
    
    A(m,l) = A(m,l) + Atemp(2,1);
    A(m,m) = A(m,m) + Atemp(2,2);
    A(m,n) = A(m,n) + Atemp(2,3);
    
    A(n,l) = A(n,l) + Atemp(3,1);
    A(n,m) = A(n,m) + Atemp(3,2);
    A(n,n) = A(n,n) + Atemp(3,3);
    
    % Compute B matrix and assemble it
    Btemp = detJ*BtempC;
    B(l,l) = B(l,l) + Btemp(1,1);
    B(l,m) = B(l,m) + Btemp(1,2);
    B(l,n) = B(l,n) + Btemp(1,3);
    
    B(m,l) = B(m,l) + Btemp(2,1);
    B(m,m) = B(m,m) + Btemp(2,2);
    B(m,n) = B(m,n) + Btemp(2,3);
    
    B(n,l) = B(n,l) + Btemp(3,1);
    B(n,m) = B(n,m) + Btemp(3,2);
    B(n,n) = B(n,n) + Btemp(3,3);
end

end