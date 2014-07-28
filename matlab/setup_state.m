%sets up the values for the global state variable
function [state,problem] = setup_state(state,problem,x0f,x0c)

    %Populate these values
    problem.n_constrained = problem.n-problem.n_free;
    %allocate the space for the working vectors
    state.temp1 = zeros(problem.n_constrained,1);
    state.temp2 = zeros(problem.n_constrained,1);
     
    %calculate the denominators for the stopping criteria
    state.rel_p_res = max(norm([problem.A,problem.b],'inf'),1);
    state.rel_d_res = max(norm(...
        [problem.A',[sparse(problem.n_free,problem.n_constrained);speye(problem.n_constrained)],problem.c],'inf'),1);
    state.rel_g_res = norm([problem.c;problem.b;1],'inf');
    
    %calculate the complexity parameter
    state.nu = problem.n_pos; %each soc cone is complexity one
    state.nu = state.nu + problem.n_exp_cones*3; %each of these of complexity 3

    state.sigma = NaN;
     
    %allocate the starting point
    state.y  = ones(problem.m,1);
    state.s  = ones(problem.n_constrained,1);
    state.xf = x0f;
    state.xc = x0c;
    state.tau = 1;
    state.kappa = 1;
    
    %calculate a feasible centered dual slack
    state.temp1  = eval_grad(problem,state.xc);
    state.s      = -(state.tau*state.kappa)/(state.nu+1)*state.temp1;
    state.temp1  = 1/(state.nu+1)*state.temp1;
    qtx          = state.s'*state.xc;
    vtx          = state.temp1'*state.xc;
    state.s      = state.s-(qtx/(vtx+1))*state.temp1; 
    
    %sanity check verify that the dual slack is feasible
    if(~eval_dual_feas(problem,state.s))
        fprintf('error, initial dual slack not feasible');
        return;
    end 
   
    %calculate the initial centrality and save in state.mu0
     dga        = state.xc'*state.s+state.kappa*state.tau;
     mua        = dga/(state.nu+1);
     state.mu0  = mua;
     state.mu   = mua;
    
    %initialize some counters
    state.c_iter = 0;
    state.m_iter = 0;
    state.b_iter = 0;
    state.c_backtrack_iter = 0;
    state.kkt_solves = 0;
    state.centering_iterations = 0;

    %allocate space for the residuals and calculate them
    state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
    state.d_res         = -problem.c*state.tau+problem.A'*state.y;
    state.d_res(problem.n_free+1:problem.n) = state.d_res(problem.n_free+1:problem.n) + state.s;

    if(~isempty(state.xf)) %if state.xf is empty then the c'x would result in an empty matrix
        cfxf = problem.c(1:problem.n_free)'*state.xf;
    else
        cfxf = 0;
    end  
    state.g_res         = - problem.b'*state.y  +cfxf+problem.c(problem.n_free+1:problem.n)'*state.xc + state.kappa; 
    
    %calculate the relative gap
    state.ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
    state.bty           = problem.b'*state.y;
    state.relative_gap  = abs( state.ctx - state.bty )/( state.tau + abs(state.bty) );
     
    %calculate the residual norms
    state.n_p_res       = norm(state.p_res,'inf')/state.rel_p_res;
    state.n_d_res       = norm(state.d_res,'inf')/state.rel_d_res;
    state.n_g_res       = abs(state.g_res)/state.rel_g_res;;

    %Failure flag
    state.fail = false;
   
end
