%Build the linear system for the atd direction
    %------------------------------------------
    
    %Calculate the hessian
    if(pars.use_nesterov_todd_centering)
        %If Nesterov Todd centering is enabled
        %Caclulate the hessians at the nt scaling points 
        H = eval_hessian_nt(problem,state.xc,state.s); 
    else
        H = eval_hessian(problem,state.xc);
    end
    
    %Shorthand 
    m = problem.m;
    n = problem.n;
    nc = problem.n_constrained;
    nf = problem.n_free;

    %Build the system matrix
    Sc  = -[sparse(nf,nc);speye(nc)];
    K5= [[sparse(m,m) , problem.A                         ,-problem.b     ,sparse(m,nc)  ,sparse(m,1)];...
         [-problem.A' , sparse(n,n)                       ,problem.c      ,Sc           ,sparse(n,1)];...
         [problem.b'  , -problem.c'                       ,sparse(1,1)    ,sparse(1,nc)  ,-1         ];...
         [sparse(nc,m), sparse(nc,nf) ,sparse(state.mu*H) ,sparse(nc,1)   ,speye(nc,nc)  ,sparse(nc,1)];...
         [sparse(1,m) , sparse(1,n)                       ,state.kappa    ,sparse(1,nc)  ,state.tau   ]];
          
          
    %Build the rhs
    rhs = [state.p_res;state.d_res;state.g_res;-state.s;-state.tau*state.kappa];  

    %Solve for the direction (The magic of UMFPACK)
    state.d       = K5\rhs;
    state.dy      = state.d(1:m);
    state.dxf     = state.d(m+1:m+nf);
    state.dxc     = state.d(m+nf+1:m+n);
    state.dtau    = state.d(m+n+1);
    state.ds      = state.d(m+n+2:m+n+nc+1);
    state.dkappa  = state.d(m+n+nc+2);

    clear 'K5' 'rhs' 'state.d'

