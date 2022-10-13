clear
clc;

global A B

% Generating a random problem:
% Minimize 1/2 * |AX-B|_F^2
%   s.t.   tr(X) = 1
%          X >= 0 

A = [];
B = [];

n = 2000;
m = 4000;
q = 10;

beta = 0.0;

rng(123456);

density = 10^(-4);

A = sprand(m,n,density);

[i,j,v] = find(A);

for k = 1:length(i)
    A(i(k),j(k)) = 2 * A(i(k),j(k)) - 1;
end

pos = randperm(n,2*q);

Z = sparse(n,n);
for i = 1:q
    num = -pi + 2.0 * pi * rand;
    u = sparse([pos(2*i-1) pos(2*i)],[1 1],[cos(num) sin(num)],n,1);
    Z = Z + u * u';
end
    
Z = ( Z + Z' ) / 2.0;

nnz(Z);

B = A * Z;

L = norm(A'*A,'fro');

% Set the starting point

X0 = ( 1.0 - beta )/n * speye(n) + beta * sparse(1,1,1.0,n,n);

% solver = 1: Gradient inexact projection method employing constant step size
% solver = 2: Gradient exact projection method employing constant step size
% solver = 3: Gradient inexact projection method employing Armijo step size
% solver = 4: Gradient exact projection method employing Armijo step size

solver = 1;

% Set the error tolerance option (see evalphi.m file)

erropt = 5;

% Call the solver

if ( solver == 1 ) 
    [X,f,time,outiter,info] = GInexPM(n,X0,L,erropt);
end

if ( solver == 2 ) 
    [X,f,time,outiter,info] = GExactPM(n,X0,L);
end

if ( solver == 3 ) 
    [X,f,time,outiter,nfev,info] = GInexPMArmijo(n,X0);
end

if ( solver == 4 ) 
    [X,f,time,outiter,nfev,info] = GExactPMArmijo(n,X0);
end