clear all
    %build a random entropy problem
    m = 10;
    n = 20;
    A = randn(m,n);
    b = A*ones(n,1);
    AA= [sparse(m+n,n),[sparse(m,n);speye(n)],[A;sparse(n,n)]];
    bb= [b;ones(n,1)];
    cc= -[ones(n,1);zeros(2*n,1)];
    N = 3*n;
    M = n+m;
    
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
    
    % starting point:
    u0  = -ones(n,1);
    v00 = ones(n,1);  
    x0  = 0.5*ones(n,1);
    
    x0c  = [u0;v00;x0];
    x0f        = []; 
    nscs 

