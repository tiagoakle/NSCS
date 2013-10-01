%NSCS matlab minimal version

clear all
clc
%Require
% in structure problem:

%m,n,p,x0,d0,
%n_free
%n_pos
%n_soc_cones
%soc_cones[]
%n_sdp_cones
%sdp_cones[]
%n_exp_cones

%Starting point
%x0f, x0c;

%Problem definition in problem strucutre
% A
% b
% c



%------------- Define a testing problem -----------
%Build an lp problem
m = 50;
n = 100;
nf = 0;
nc = n-nf;

A  = randn(m,n);
x  = [randn(nf,1);rand(nc,1)];
b  = A*x;

c          = A'*randn(m,1);
c(nf+1:n)  = c(nf+1:n) + rand(n-nf,1);

%Generate a random start
x0 = [randn(nf,1);rand(nc,1)];
problem = struct;
problem.A = A;
problem.b = b;
problem.c = c;
problem.m = m;
problem.n = n;
problem.n_free = nf;
problem.n_pos = nc;
problem.soc_cones = 0;
problem.n_soc_cones = 0;
problem.n_sdp_cones = 0;
problem.sdp_cones     = 0;
problem.n_exp_cones   = 0;
problem.n_power_cones = 0;
x0f        = randn(nf,1);
x0c        = rand(nc,1);

clear 'm' 'n' 'nf' 'nc' 'A' 'b' 'c' 'x';


%--------- End of problem definition -------------

%Set up the default parameters
    pars.max_iter   = 30;
    pars.max_affine_backtrack_iter = 100;
    pars.max_centering_backtrack_iter = 100;
    pars.max_c_iter = 20;
    pars.backtrack_affine_constant = 0.5;
    pars.backtrack_centering_constant = 0.5;
    pars.beta       = 0.9;
    pars.theta      = 0.8;
    pars.large_neighborhood = 1.001;
    pars.eta        = 0.995;
    pars.use_nesterov_todd_centering = false;
    pars.stop_primal= 1e-10;
    pars.stop_dual  = 1e-10;
    pars.stop_gap   = 1e-10;

%initialize the state
state = struct;
state.xc = x0c;
state.xf = x0f;

%------------------------------------------------
% Validate problem input and initialize the state
%------------------------------------------------

%Make sure the number of constraints adds to p;
tot_constraints = problem.n_free+problem.n_pos+sum(problem.soc_cones)+sum(problem.sdp_cones)+...
                  problem.n_exp_cones*3+problem.n_power_cones*3;

if(tot_constraints ~= problem.n)
    fprintf('Error, sum of all constraints must equal n\n');
    return;
end

%Populate these values
problem.n_constrained = problem.n-problem.n_free;

%Verify that the inital d is feasible
if(~eval_primal_feas(problem,state.xc))
    fprintf('Error, initial primal slack not feasible');
    return;
end

%----------------------------------------------
% Print the header
%-----------------------------------------------
fprintf('%2s  %2s  %6s   %6s     %6s     %6s       %6s     %6s     %6s\n',...
                         'it','cit',...
                         'a_a','mu',...
                         'tau','kap',...
                         'p_res','d_res','g_res');



%-----------------------------------------------
%Initialize the variables that are a function of 
% the initial point and parameters
%-----------------------------------------------

%Calculate the denominators for the stopping criteria
state.rel_p_res = max(norm([problem.A,problem.b],'inf'),1);
state.rel_d_res = max(norm(...
    [problem.A',[sparse(problem.n_free,problem.n_constrained);speye(problem.n_constrained)],problem.c],'inf'),1);
state.rel_g_res = norm([problem.c;problem.b;1],'inf');

%Calculate the complexity parameter
state.nu = problem.n_pos+problem.n_soc_cones*2;
state.nu = state.nu + sum(problem.sdp_cones);
state.nu = state.nu + problem.n_exp_cones*3+problem.n_power_cones*3;

% Complete the starting point
state.y  = ones(problem.m,1);
state.s  = ones(problem.n_constrained,1);
state.xf = x0f;
state.xc = x0c;
state.tau = 1;
state.kappa = 1;

%Allocate the space for the working vectors
state.temp1 = zeros(problem.n_constrained,1);
state.temp2 = zeros(problem.n_constrained,1);

%Caclculate a feasible centered dual slack
state.temp1  = eval_grad(problem,state.xc);
state.s      = -(state.tau*state.kappa)/(state.nu+1)*state.temp1;
state.temp1  = 1/(state.nu+1)*state.temp1;
qtx          = state.s'*state.xc;
vtx          = state.temp1'*state.xc;
state.s      = state.s-(qtx/(vtx+1))*state.temp1;

%Calculate the initial centrality
 dga        = state.xc'*state.s+state.kappa*state.tau;
 mua        = dga/(state.nu+1);

state.mu    = mua;

%initialize some counters
state.c_iter = 0;
state.m_iter = 0;
state.b_iter = 0;
state.c_backtrack_iter = 0;

%Shorthand symbols
n           = problem.n;
nc          = problem.n_constrained;
nf          = problem.n_free;
m           = problem.m;

%Allocate space for the residuals
state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
state.d_res         = -problem.c*state.tau+problem.A'*state.y;
state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
state.g_res         = - problem.b'*state.y  + problem.c(1:nf)'*state.xf+problem.c(nf+1:n)'*state.xc + state.kappa; 
%Calculate the residual norms
state.n_p_res       = norm(state.p_res,'inf')/state.rel_p_res;
state.n_d_res       = norm(state.p_res,'inf')/state.rel_d_res;
state.n_g_res       = abs(state.g_res)/state.rel_g_res;;

%Print the iteration log
%iter centering iter, 
fprintf('%2i  %2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e\n',...
                     state.m_iter,state.c_iter,...
                     0,state.mu,...
                     state.tau,state.kappa,...
                     state.n_p_res,state.n_d_res,state.n_g_res);


%-----------------------------------------------
%Start of the major iteration
%-----------------------------------------------
m_iter = 0;
for m_iter = 0:pars.max_iter
   
    state.m_iter = m_iter;
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
    %  Start of the approximate tangent direction linesearch
    %----------------------------------------------------------
    state.a_affine  = 1.0;
    %Find the largest step to the boundary of tau, kappa 
    % backtrack until the centering measure is smaller than mu*pars.theta

    if(state.dtau<0)
        state.a_affine = min(state.a_affine,-state.tau/state.dtau);
    end
    if(state.dkappa<0)
        state.a_affine = min(state.a_affine,-state.kappa/state.dkappa);
    end
    %Take a multiple of the maximum length
    state.a_affine = state.a_affine*pars.eta;
    
    % Backtrack loop
    b_iter = 0;
    for b_iter = 1:pars.max_affine_backtrack_iter
        state.b_iter = b_iter;

        %Evaluate the trial point
        ya         = state.y         + state.a_affine*state.dy; 
        xaf        = state.xf        + state.a_affine*state.dxf;
        xac        = state.xc        + state.a_affine*state.dxc;
        taua       = state.tau       + state.a_affine*state.dtau;
        sa         = state.s         + state.a_affine*state.ds;
        kappaa     = state.kappa     + state.a_affine*state.dkappa;

        %Check if the present point is primal dual feasible
        p_feas     = eval_primal_feas(problem,xac);
        if(p_feas)
            d_feas     = eval_dual_feas(problem,sa);
            if(d_feas)
                %Calculate mu and the duality gap at the trial point
                dga        = xac'*sa+kappaa*taua;
                mua        = dga/(state.nu+1);
                %Evaluate the centrality measure ||sa+mua*g(xa)||_H^{-1}(xa)
                H           = eval_hessian(problem,xac);
                state.g     = eval_grad(problem,xac); %we can save this with a product..
                state.temp1 = +mua*state.g;
                state.temp1 = sa+state.temp1;
                state.temp2 = H\state.temp1;
                state.cent  = sqrt(state.temp1'*state.temp2);

                %XXX: DEBUG
                fprintf('B_iter %i, Cent %g large_neigh*mu*theta: %g\n',b_iter,state.cent,pars.large_neighborhood*mua*pars.theta);
                
                if(state.cent <= pars.large_neighborhood*mua*pars.theta) %Stop when this is satisfied
                    break;
                end
            else
                %not dual infeasible
            end
        else
            %not primal infeasible 
        end 
        state.a_affine  = state.a_affine*pars.backtrack_affine_constant;
    end %End of backtrack loop

    %If the maximum number of iterates was reached report an error and exit
    if(state.cent > pars.large_neighborhood*mua*pars.theta)
        fprintf('Backtracking line search failed is backtrack_affine_constant too large?\n');
        state.exit_reason = 'affine backtrack line search fail';
        break;
    end

    %If the line-search succeeded take the step
    state.y     = ya;
    state.xf    = xaf;
    state.xc    = xac;
    state.tau   = taua;
    state.s     = sa;
    state.kappa = kappaa;
    state.mu    = mua;
   
    clear 'ya' 'xaf' 'xac' 'taua' 'sa' 'kappaa' 'mua'
    %          Centering process
    %----------------------------------------------------------
    fprintf('Cent before centering: %g mu: %g large_neigh*mu*theta %g\n',state.cent,state.mu,...
    state.mu*pars.theta*pars.large_neighborhood);    

     for c_iter = 1:pars.max_c_iter
       state.c_iter = c_iter;
        %Evaluate the gradient and hessian
        state.g = eval_grad(problem,state.xc);
        H       = eval_hessian(problem,state.xc);

        %Build the system matrix
        Sc  = -[sparse(nf,nc);speye(nc)];
        K5= [[sparse(m,m) , problem.A                         ,-problem.b     ,sparse(m,nc)  ,sparse(m,1)];...
             [-problem.A' , sparse(n,n)                       ,problem.c      ,Sc           ,sparse(n,1)];...
             [problem.b'  , -problem.c'                       ,sparse(1,1)    ,sparse(1,nc)  ,-1         ];...
             [sparse(nc,m), sparse(nc,nf) ,sparse(state.mu*H) ,sparse(nc,1)   ,-speye(nc,nc)  ,sparse(nc,1)];...
             [sparse(1,m) , sparse(1,n)                       ,state.kappa    ,sparse(1,nc)  ,state.tau   ]];
        
        
        %Evaluate the centering condition
        state.temp1 = state.s+state.mu*state.g;
%XXX DEBUG
        state.cent2norm = norm(state.temp1);
%END DEBUG
        state.temp2 = H\state.temp1;
        state.cent  = sqrt(state.temp1'*state.temp2);
      
   
        %Build the rhs
        rhs                 = zeros(n+m+nc+2,1);
        rhs(m+n+2:m+n+1+nc) = -state.s-state.mu*state.g;
        rhs(m+n+nc+2)       = state.mu - state.tau*state.kappa;

        %Solve for the direction 
        state.d    = K5\rhs;
        state.dy   = state.d(1:m);
        state.dxf  = state.d(m+1:m+nf);
        state.dxc  = state.d(m+nf+1:m+n);
        state.dtau = state.d(m+n+1);
        state.ds   = state.d(m+n+2:m+n+nc+1);
        state.dk   = state.d(m+n+nc+2);
        
        %Centering backtrack 
        state.a_cent  = 1.0;
        %Find the largest step to the boundary of tau, kappa 
        % backtrack until the centering measure is smaller than mu*pars.theta
        
        if(state.dtau<0)
            state.a_cent = min(state.a_cent,-state.tau/state.dtau);
        end
        if(state.dkappa<0)
            state.a_cent = min(state.a_cent,-state.kappa/state.dkappa);
        end
        
        %Take a multiple of the maximum length
        state.a_cent = state.a_cent*pars.eta;
        %state.cent = Inf;
        %Backtracking iteration

        %fprintf('Centrality at bck: %g\n',state.cent);    
        c_backtrack_iter = 0;
        for c_backtrack_iter = 1:pars.max_centering_backtrack_iter
            state.c_backtrack_iter = c_backtrack_iter;

            %Calculate the trial point
            xca        = state.xc        + state.a_cent*state.dxc;
            taua       = state.tau       + state.a_cent*state.dtau;
            sa         = state.s         + state.a_cent*state.ds;
            kappaa     = state.kappa     + state.a_cent*state.dkappa;
            
            %Calculate mu and the duality gap at the trial point
            dga        = xca'*sa+kappaa*taua;
            mua        = dga/(state.nu+1);
           
            %Evaluate the centrality measure ||sa+mua*g(xa)||_H^{-1}(xa)
            H           = eval_hessian(problem,xca);
            state.g     = eval_grad(problem,xca); %we can save this with a product..
            state.temp1 = sa+mua*state.g;
            state.temp2 = H\state.temp1;
            cent        = sqrt(state.temp1'*state.temp2);
            
            %XXX DEBUG
            cent2norm   = norm(state.temp1);
            if(m_iter == 1)
                if(~exist('logs','var'))
                    logs = zeros(pars.max_centering_backtrack_iter,3);
                end  
                    logs(c_backtrack_iter,1) = state.a_cent;
                    logs(c_backtrack_iter,2) = cent;
                    logs(c_backtrack_iter,3) = cent2norm;
            end
            %END DEBUG 
            if(cent < state.cent) %Accept the step if it decreases the centrality... 
            %XXX: add an armijo criteria to this??
                state.cent = cent;
                break;
            else
                state.a_cent = state.a_cent*pars.backtrack_centering_constant;
            end
                     
        end %End of backtracking iteration

        state.y     = state.y         + state.a_cent*state.dy;
        state.xf    = state.xf        + state.a_cent*state.dxf;
        state.xc    = state.xc        + state.a_cent*state.dxc; 
        state.tau   = state.tau       + state.a_cent*state.dtau;
        state.s     = state.s         + state.a_cent*state.ds;
        state.kappa = state.kappa     + state.a_cent*state.dkappa;
        state.mu    = mua;
        
        %fprintf('Cent iter :%i bkt %i, centrality %g \n',c_iter, state.c_backtrack_iter,state.cent)

        if(state.cent <= mua*pars.theta) %Stop when the centering measure is smaller than beta
            break;
        end

  end %End of centering iteration
 
  %Check if the centering iteration failed   
  if(state.cent > state.mu*pars.theta)
        fprintf('Centering failed \n');
        state.exit_reason = 'Centering fail';
        break; 
  end

  fprintf('Cent after centering: %g mu: %g large_neigh*mu*theta %g\n',state.cent,state.mu,...
    state.mu*pars.theta*pars.large_neighborhood);    

    %Calculate the residuals
    state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
    state.d_res         = -problem.c*state.tau+problem.A'*state.y;
    state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
    state.g_res         = - problem.b'*state.y +...
                            problem.c(1:nf)'*state.xf+problem.c(nf+1:n)'*state.xc + state.kappa; 
    
    %Calculate the residual norms
    state.n_p_res       = norm(state.p_res,'inf')/state.rel_p_res;
    state.n_d_res       = norm(state.p_res,'inf')/state.rel_d_res;
    state.n_g_res       = abs(state.g_res)/state.rel_g_res;;

    %Print the iteration log
    %iter centering iter, 
    fprintf('%2i  %2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e\n',...
                         state.m_iter,state.c_iter,...
                         state.a_affine,state.mu,...
                         state.tau,state.kappa,...
                         state.n_p_res,state.n_d_res,state.n_g_res);

    %Evaluate the stopping criteria
    if(state.n_p_res < pars.stop_primal && state.n_p_res < pars.stop_dual && state.n_g_res < pars.stop_gap)
        if(state.kappa<state.tau-eps)
            state.exit_reason  = 'Optimal';
        elseif (state.kappa>state.ta+eps)
            state.exit_reason  = 'Infeasible';
            %XXX detect p.d inf
        else
            state.exit_reason = 'Ill Posed';
        end
        break;
    end
 

end %End of main loop


