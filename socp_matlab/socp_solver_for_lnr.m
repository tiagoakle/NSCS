%Exercise, an SOCP solver for least norm of x subject to measurement noise
% b = Ax + n where n is gausian noise
% minimize ||x||^2 st ||Ax-b||<sigma^2

%This is an SOCP of the following form

%minimze x_0 
% subject to 
% [0,A,0,I][x_0]  = [b      ]
% [0,0,1,0][x  ]    [sigma^2]
%          [r_0]
%          [r  ]

% [x_0,x] \in L_{n+1} [r_0,r] \in L_{n+1}

%----------------------------------------
%Problem parameters
m = 10;
n = 15;
sigma = 0.1;

%Build a random problem
A = randn(m,n);
x_exact = randn(n,1);
b = A*x_exact+1/sqrt(m)*sigma*randn(m,1);

problem = struct;
problem.soc_cones = zeros(2,1);
problem.n_soc_cones = 2;
problem.soc_cones(1) = n+1;
problem.soc_cones(2) = n+1;

%Find an initial point 
x = ones(n,1);
x_0 = norm(x) + 1;
r = ones(n,1);
r_0 = norm(r) + 1;


