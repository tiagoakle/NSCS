clear all
addpath ../../coneopt
%Create a random entropy problem and
% solve it with coneopt
    
   m = 100;
   n = 200;
   A = sprandn(m,n,0.5);
   b = A*ones(n,1);
   AA= [sparse(m+n,n),[sparse(m,n);speye(n)],[A;sparse(n,n)]];
   bb= [b;ones(n,1)];
   cc= -[ones(n,1);zeros(2*n,1)];
   N = 3*n;
   M = n+m;
 
 % build cone:
    K.npos = 0;
    K.npow = 0;
    K.nexp = n;
    K.nlog = 0;
    K      = getbarrpar(K);

    pars.m    = M;
    pars.n    = N;
    pars.centmeastype = 5; 
    pars.echo = 4;
    pars.secord = 0;


  % starting point:
  u0  = -ones(n,1);
  v00 = ones(n,1);  
  x0  = 0.5*ones(n,1);
  
  v0.x  = [u0;v00;x0];

% call to coneopt:
tic
R = coneopt(AA,bb,cc,v0,K,pars);
coneopt_time = toc

% call nscs
%Extract the problem data and build the problem structure
    problem = struct;
    problem.A = AA;
    problem.b = bb;
    problem.c = cc;
    %Problem parameters
    problem.m = M;
    problem.n = N;
    problem.n_free = 0;
    problem.n_constrained = N;
    problem.n_pos       = 0;
    problem.soc_cones   = 0;
    problem.n_soc_cones = 0;
    problem.n_sdp_cones = 0;
    problem.sdp_cones     = 0;
    problem.n_exp_cones   = n;
    problem.n_power_cones = 0;
   
x0c = [u0;v00;x0];
x0f = [];
tic 
nscs
nscs_time = toc;

%Now call sedumi via CVX
% Entropy maximization
tic 
cvx_begin
    cvx_solver('sedumi')
    variable x(n/3)
    maximize sum(entr(x))
    A*x == b
cvx_end
sedumi_time = toc;


