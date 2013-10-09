clear all

    addpath '../'
    %Load the file that contains the indices for the
    %ufget netlib lps which are in standard form
    load 'standard_form_indices.mat' 

    %Choose a problem from the list
    problem_index = 14; 
    %Extract the problem 
    problem_uf_ix = st_ix(problem_index);
    %Get the problem from ufget
    P = UFget(problem_uf_ix);
      
    %Extract the name
    prob_name = [P.name];
    %Substitute front slash for space
    prob_name(find(prob_name=='/'))=' ';
    
    %Extract the problem data and build the problem structure
    problem = struct;
    problem.A = P.A;
    problem.b = P.b;
    problem.c = P.aux.c;
    %Problem parameters
    problem.m = size(problem.A,1);
    problem.n = size(problem.A,2);
    problem.n_free = 0;
    problem.n_constrained = problem.n;
    problem.n_pos       = problem.n;
    problem.soc_cones   = 0;
    problem.n_soc_cones = 0;
    problem.n_sdp_cones = 0;
    problem.sdp_cones     = 0;
    problem.n_exp_cones   = 0;
    problem.n_power_cones = 0;

    pars = struct;
    pars.delta = 1.e-10 ;
    pars.gamma = 1.e-10 ;
    pars.max_iter_ref_rounds = 3;

    x0c        = 1*rand(problem.n,1);
    fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,problem.m,problem.n);

    %Generate a dual value for x
    s        = eval_grad(problem,x0c);
    mu = 1;
    kappa = 1;
    tau  = 1;
    %Evaluate the hesian at the point
    H = eval_hessian(problem,x0c);
   
   %Build the rhs
   r1 = zeros(problem.m,1);
   r2 = zeros(problem.n,1);
   r3 = 0;
   r4 = 0;  %Remember coneopt reverses these
   r5 = -s;
   
    %Call the solver 
    d_mixed = solve_linear_system_mixed(H,mu,kappa,tau,problem,pars,r1,r2,r3,r4,r5);
   
    %----------------------------------------------------------------------
    % Compare to coneopt's solver
    %----------------------------------------------------------------------
 
    m = problem.m;
    n = problem.n;
    A = problem.A;
    b = sparse(problem.b);
    c = sparse(problem.c);
   
    %Make a struct for the parameters that slvhomkkt uses
    parssolve = struct;
    parssolve.permuteM = false;
    parssolve.cholinc = true;
    
    [d_coneopt,CF] = slvhomkkt(H,mu,A,b,c,tau,kappa,...
    r1,r2,r3,r4,r5,parssolve);
    state.dx     = d_coneopt{1};
    state.dtau   = d_coneopt{2};
    state.dy     = d_coneopt{3};
    state.ds     = d_coneopt{4};
    state.dkappa = d_coneopt{5};   

    %print the difference in norm
    fprintf('ndx %g, ndy %g, nds %g ,ndt %g, ndk %g \n',norm(state.dx-d_mixed.dxc),...
                                                             norm(state.dxc-d_mixed.dy),...
                                                             norm(state.ds-d_mixed.ds),...
                                                             abs(state.tau-d_mixed.tau),...
                                                             abs(state.kappa-d_mixed.kappa));


