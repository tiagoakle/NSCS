function [x] = minentropy_nscs(A,b,d,varargin)

% min_x   sum_j  d_j x_j log(x_j) 
% s.t.     Ax = b
%           x >= 0
%

[M,N]   = size(A);
b       = b(:);
M2      = length(b);

if M ~= M2
    error('A and b do not match');
end

problem.n = 3*N;
problem.m = N+M;

problem.n_free = 0;
problem.n_constrained = 3*N;
problem.n_pos = 0;
problem.soc_cones   = 0;
problem.n_soc_cones = 0;
problem.n_exp_cones   = N;

% build A:
nnzA    = nnz(A);
AA      = sparse([],[],[],problem.m,problem.n,nnzA+N);
tmp     = N*(N+M);
tmp     = tmp + 1:(N+M+1):(2*N*(N+M)-M);
AA(tmp) = 1;
tmp     = N+1:problem.m;
tmp2    = 2*N+1:problem.n;
AA(tmp,tmp2) = A;

% build b:
bb             = ones(problem.m,1);
bb(N+1:problem.m) = b;
bb             = full(bb);

% build c:
cc      = zeros(problem.n,1);
cc(1:N) = -d;
cc      = sparse(cc);

problem.A = AA;
problem.b = bb;
problem.c = cc;

% starting point:
u0  = -ones(N,1);
v00 = ones(N,1);  
x0  = 0.5*ones(N,1);

x0c = [u0;v00;x0];
x0f = [];
pars = set_default_pars_nscs_long_step(); 
pars.second_order = false;
pars.neigh = 1.0;

pars.stop_primal= 1e-6;                 %stopping criteria p_res/rel_p_res<stop_primal.
pars.stop_dual  = 1e-6;
pars.stop_gap   = 1e-6;
pars.stop_mu    = 1e-8;
pars.stop_tau_kappa = 1.e-5;
pars

 
[xc,xf,y,s,t,k,info] = nscs_long_step(problem,x0f,x0c,pars);
x = xc;

end
