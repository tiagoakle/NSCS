%Test dual feasibility LP

%Problem parameters
problem.n = 10;
problem.n_free = 0;
problem.n_constrained = 10;
problem.n_pos       = 10;
problem.soc_cones   = 0;
problem.n_soc_cones = 0;
problem.n_sdp_cones = 0;
problem.sdp_cones     = 0;
problem.n_exp_cones   = 0;
problem.n_power_cones = 0;


x = rand(10,1);
if(~eval_dual_feas(problem,x)) error('Test dual feas LP failed'); end
x(1) = 0;
if(eval_dual_feas(problem,x)) error('Test dual feas LP failed feas = true and x(1) = 0'); end

%Test dual feas exp
problem.n_exp_cones = 10;
problem.n_pos       = 0;
w00 = -ones(10,1);
v00 = ones(10,1);  
y00 = 0.5*ones(10,1);
x   = [w00,v00,y00];
if(~eval_primal_feas(problem,x)) error('Test primal feas exp failed'); end
x   = [-y00,v00/exp(1),-w00];
if(~eval_dual_feas(problem,x)) error('Test dual feas exp failed'); end



disp('No errors');
