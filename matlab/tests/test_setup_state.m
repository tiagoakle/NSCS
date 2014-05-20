msg = 'Test setup state: ';

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
 
state = struct;
[state,prob] = setup_state(state,problem,x0f,x0c);

msg2=' computing number of constrained vars';
assert(problem.n_constrained  == problem.n,[msg,msg2]);

msg2=' size of allocated work vectors';
assert(length(state.temp1)==problem.n_constrained,[msg,msg2]);
assert(length(state.temp2)==problem.n_constrained,[msg,msg2]);

msg2='Initial denominators';
assert(state.rel_p_res == 2,[msg,msg2]);
assert(state.rel_d_res == 3,[msg,msg2]);
assert(state.rel_g_res == 1,[msg,msg2]);

msg2='Initial complexity';
assert(state.nu == 40,[msg,msg2]);

msg2='Initial point';
assert(all(state.y == ones(problem.m,1)),[msg,msg2]); 

%check that the s is centered
gx0 = eval_grad(problem,x0c);
assert(norm([state.s;state.tau]+[state.mu*gx0;-state.mu/state.kappa])<1.e-14,[msg,msg2]);

assert(all(state.xc == x0c),[msg,msg2]); 
assert(all(state.xf == x0f),[msg,msg2]); 
assert(state.tau   ==1,[msg,msg2]); 
assert(state.kappa ==1,[msg,msg2]); 


msg2= 'Initial centrality ';
assert(state.mu==(state.xc'*state.s+state.kappa*state.tau)/(state.nu+1),[msg,msg2]);

msg2='Set counters to zero';
assert(state.c_iter == 0,[msg,msg2])
assert(state.m_iter == 0,[msg,msg2])
assert(state.b_iter == 0,[msg,msg2])
assert(state.c_backtrack_iter == 0,[msg,msg2])
assert(state.kkt_solves == 0,[msg,msg2])
assert(state.centering_iterations == 0,[msg,msg2])

msg2='Initial residuals';
assert(all(state.p_res ==  problem.b*state.tau-problem.A*[state.xf;state.xc]),[msg,msg2]);
err_d_res = -problem.c*state.tau+problem.A'*state.y;
err_d_res(problem.n_free+1:problem.n) = err_d_res(problem.n_free+1:problem.n)+ state.s;
err_d_res = err_d_res - state.d_res;
assert(norm(err_d_res)==0,[msg,msg2]);

if(problem.n_free ==0 )
    cfxf=0;
else
    cfxf= problem.c(1:problem.n_free)'*state.xf;
end

assert(state.g_res == - problem.b'*state.y  +cfxf+problem.c(problem.n_free+1:problem.n)'*state.xc + state.kappa,[msg,msg2]);
    
%calculate the relative gap
    ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
    bty           = problem.b'*state.y;
assert(state.relative_gap  == abs( ctx - bty )/( state.tau + abs(bty) ),[msg,msg2]);
     
    %calculate the residual norms
assert(state.n_p_res == norm(state.p_res,'inf')/state.rel_p_res);
assert(state.n_d_res == norm(state.d_res,'inf')/state.rel_d_res);
assert(state.n_g_res == abs(state.g_res)/state.rel_g_res);
%Check the flag is correctly set 
  assert(state.fail == false);
