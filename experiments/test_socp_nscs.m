clear all
close all
addpath ../matlab
%Test the socp vertion of nscs

n = 100; %Size of the cone
m = 50; %Number of costraints
nu = 1; %Complexity of the barrier

%-------
%Generate a feasible primal point
x = randn(n,1);
x(1) = norm(x(2:end))+rand(1);
x_feas = x;

%Generate a feasible dual point
s = randn(n,1);
s(1) = norm(s(2:end))+rand(1);
s_feas = s;

%Generate a feasible and bounded primal dual problem
A = sprandn(m,n,0.5);
b = A*x_feas;
s = rand(n,1);
c = A'*ones(m,1) + s_feas;
fprintf('Generated a random feasible SOCP with %i constraints and %i variables\n',m,n);
clear x_feas s_feas

%Choose an initial point
x0c = [1;zeros(n-1,1)];
x0f  = [];

%Extract the problem data and build the problem structure
problem = struct;
problem.A = A;
problem.b = b;
problem.c = c;

%Problem parameters
problem.m = m;
problem.n = n;
problem.n_free = 0;
problem.n_constrained = n;
problem.n_pos       = 0;
problem.n_soc_cones = 1;
problem.soc_cones   = [n];
problem.n_sdp_cones = 0;
problem.sdp_cones     = [];
problem.n_exp_cones   = 0;

%Populate the default structure 
set_default_pars_nscs_long_step;

%--------------------------------------------------------------------------
% Solve with nscs long step with nt scaling
%--------------------------------------------------------------------------
pars.use_nesterov_todd_scaling = true;
[xc,xf,y,z,info] = nscs_long_step(problem,x0f,x0c,pars);
nscs_lsnt_kkt = info.kkt_solves; 
nscs_lsnt_sta = info.exit_reason;
 
