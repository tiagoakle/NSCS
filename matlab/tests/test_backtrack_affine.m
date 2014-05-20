clear all
%Test the affine backtracking 
%Must terminate in a primal dual feasible point with enough centrality
%Must generate its iterates along state.dxc state.dtau,... etc

%Feasible primal dual problem and present iterate

%Parameters for the test
pars.max_affine_backtrack_iter = 100;
pars.neigh = 1.0;
pars.print = 4;
pars.backtrack_affine_constant = 0.5;

%Test the allocation
problem.n = 40;
problem.m = 10;
problem.A = [speye(problem.m),sparse(problem.m,problem.n-problem.m)];
problem.b = ones(problem.m,1);
problem.c = ones(problem.n,1);

problem.n = 40;
problem.n_free = 0;
problem.n_constrained = 40;
problem.n_pos       = 10;
problem.soc_cones   = 0;
problem.n_soc_cones = 0;
problem.n_sdp_cones = 0;
problem.sdp_cones     = 0;
problem.n_exp_cones   = 10;
problem.n_power_cones = 0;

x0f = [];
x0c = ones(problem.n_constrained,1);
x0c(problem.n_pos+1:problem.n_pos+problem.n_exp_cones) = -0.827838387734678*ones(problem.n_exp_cones,1);
x0c(problem.n_pos+problem.n_exp_cones+1:problem.n_pos+2*problem.n_exp_cones) = 1.290927717327200*ones(problem.n_exp_cones,1);
x0c(problem.n_pos+2*problem.n_exp_cones+1:problem.n_pos+3*problem.n_exp_cones) = 0.805102001257750*ones(problem.n_exp_cones,1);

state.xc = x0c;
state.s  = +x0c;
state.y = 0;
state.tau = 1;
state.kappa = 1;
state.a_affine = 1;
state.nu = 40+30;


%Make sure the primal linear becomes feasible
state.dxc(1:problem.n_pos,1) = -10*ones(problem.n_pos,1);
state.dxc(problem.n_pos+1:problem.n_pos+3*problem.n_exp_cones,1) = -zeros(3*problem.n_exp_cones,1);
state.dtau = -10;
state.dkappa = 0;
state.dy = zeros(problem.m,1);
state.ds  = zeros(problem.n_constrained,1);

%Santity 
dga   = state.xc'*state.s+state.kappa*state.tau;
mua   = dga/(state.nu+1);
    
eval_small_neigh(problem,state.xc,state.s,mua)
state = backtrack_affine(state,pars,problem);

assert(state.fail == false,' back into primal feas');
