%NSCS matlab minimal version
%This path is used to call the linear solver in coneopt
addpath '../../coneopt/'
%clear all
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

%Quiet the matlab warning about bad scaling
warning('off','MATLAB:nearlySingularMatrix');


%Set up the default parameters
    pars.max_iter   = 100;  %Maximum outer iterations
    pars.max_affine_backtrack_iter = 300;    %Maximum affine backtracking steps
    pars.max_centering_backtrack_iter = 300; %Maximum centering backtracking steps
    pars.max_c_iter = 50;                    %Maximum centering iterations per affine  iteration
    pars.backtrack_affine_constant = 0.94;   %Affine backtracking constant
    pars.backtrack_centering_constant = 0.5; %Centering backtracking constant
    pars.beta       = 0.2;                   %Stop centering when ||dx||_H<beta
    pars.theta      = 0.8;                   %Take an affine stem if ||Psi^+||<theta*mu+
    pars.eta        = 0.9995;                %Multiple of step to the boundary
    pars.use_nesterov_todd_centering = false; %Use centering points for symmetric cones
    pars.stop_primal= 1e-6;                 %Stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual  = 1e-6;
    pars.stop_gap   = 1e-6;
    pars.stop_mu    = 1e-6;
    pars.stop_tau_kappa = 1.e-6;
    pars.solve_second_order = true;

    pars.print      = 1;                     %Level of verbosity from 0 to 11
    %Regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 100;
    pars.linear_solver = 'mixed';
    pars.centrality_measure = 1;

%------------------------------------------------
% Validate problem input and initialize the state
%------------------------------------------------

state = struct;
state.xc = x0c;
state.xf = x0f;

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
if(pars.print >0)
    fprintf(' Problem size (%i,%i) nnz(A): %i \n',problem.m,problem.n,nnz(problem.A));
    fprintf(' Free: %i, Positive %i, SOCP cones %i, Matrix %i, Exponential Cones %i\n',...
        problem.n_free,problem.n_pos,problem.n_soc_cones,...
        problem.n_sdp_cones,problem.n_exp_cones);
    fprintf(' Linear solver %s centrality %i \n',pars.linear_solver,pars.centrality_measure);
    fprintf('==========================================================================\n');
    fprintf('%2s  %2s  %6s   %6s     %6s     %6s       %6s       %6s     %6s     %6s\n',...
                         'it','cit',...
                         'a_a','mu',...
                         'tau','kap',...
                         'p_res','d_res','rel_gap','g_res');
    fprintf('==========================================================================\n');
end



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

%Allocate the space for the working vectors
state.temp1 = zeros(problem.n_constrained,1);
state.temp2 = zeros(problem.n_constrained,1);

% Complete the starting point
state.y  = ones(problem.m,1);
state.s  = ones(problem.n_constrained,1);
state.xf = x0f;
state.xc = x0c;
state.tau = 1;
state.kappa = 1;

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
 state.mu0  = mua;
 state.mu   = mua;

%initialize some counters
state.c_iter = 0;
state.m_iter = 0;
state.b_iter = 0;
state.c_backtrack_iter = 0;
state.kkt_solves = 0;
state.centering_iterations = 0;

%Shorthand symbols
n           = problem.n;
nc          = problem.n_constrained;
nf          = problem.n_free;
m           = problem.m;

%Allocate space for the residuals and calculate them
state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
state.d_res         = -problem.c*state.tau+problem.A'*state.y;
state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
if(~isempty(state.xf)) %If state.xf is empty then the c'x would result in an empty matrix
    cfxf = problem.c(1:nf)'*state.xf;
else
    cfxf = 0;
end  
state.g_res         = - problem.b'*state.y  +cfxf+problem.c(nf+1:n)'*state.xc + state.kappa; 

state.ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
state.bty           = problem.b'*state.y;
state.relative_gap  = abs( state.ctx - state.bty )/( state.tau + abs(state.bty) );


%Calculate the residual norms
state.n_p_res       = norm(state.p_res,'inf')/state.rel_p_res;
state.n_d_res       = norm(state.d_res,'inf')/state.rel_d_res;
state.n_g_res       = abs(state.g_res)/state.rel_g_res;;

%Print the iteration log
if(pars.print > 0)
    fprintf('%2i  %2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e\n',...
                     state.m_iter,state.c_iter,...
                     0,state.mu,...
                     state.tau,state.kappa,...
                     state.n_p_res,state.n_d_res,state.relative_gap,state.n_g_res);
end


%-----------------------------------------------
%Start of the major iteration
%-----------------------------------------------
m_iter = 0;
state.exit_reason = 'Max Iter Reached';
for m_iter = 1:pars.max_iter
   
    state.m_iter = m_iter;
    
    %-------------------------
    %Solve the linear system
    %-------------------------
    %Evaluate the hessian either at X or at the centering point
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
      
    %Solve the linear system 
    if(strcmp(pars.linear_solver,'umfpack')==1)
        %Assemble the matrix
        Sc  = [sparse(nf,nc);-speye(nc)];
        K5= [[sparse(m,m) , problem.A                           ,-problem.b     ,sparse(m,nc)  ,sparse(m,1)];...
             [-problem.A' , sparse(n,n)                         ,problem.c      ,Sc            ,sparse(n,1)];...
             [problem.b'  , -problem.c'                         ,sparse(1,1)    ,sparse(1,nc)  ,-1         ];...
             [sparse(nc,m), [sparse(nc,nf) ,sparse(state.mu*H)] ,sparse(nc,1)   ,speye(nc,nc)  ,sparse(nc,1)];...
             [sparse(1,m) , sparse(1,n)                         ,state.kappa    ,sparse(1,nc)  ,state.tau   ]];
         
              
        %Build the rhs
        rhs = [state.p_res;state.d_res;state.g_res;-state.s;-state.tau*state.kappa];  

        d       = K5\rhs; 
        state.dy       = d(1:m);
        state.dxf      = d(m+1:m+1+nf);
        state.dxc      = d(m+nf+1:m+n);
        state.dtau     = d(m+n+1);
        state.ds       = d(m+n+2:m+n+1+nc);
        state.dkappa   = d(m+n+nc+2);

    elseif(strcmp(pars.linear_solver,'mixed')==1)

        [d, factorization]  = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
                             state.p_res,state.d_res,state.g_res,-state.tau*state.kappa,-state.s,[]);
 
        state.dy       = d.dy;
        state.dxf      = d.dxf;
        state.dxc      = d.dxc;
        state.dtau     = d.dtau;
        state.ds       = d.ds;
        state.dkappa   = d.dkappa;
        
        %Evaluate the second order direction
        %And eventually the Mehrota predictor corrector not implemented yet
        if(pars.solve_second_order)
             
             state.g        = eval_grad(problem,state.xc);
            %d              = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
            %                 zeros(size(state.p_res)),zeros(size(state.d_res)),0,-state.dtau*state.dkappa,2*state.mu*(state.xc.^(-3)).*(state.dxc.^2)+2*(-state.s-state.ds),factorization);
            rhs             = eval_tensor(problem,state);
            d               = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
                              zeros(size(state.p_res)),zeros(size(state.d_res)),0,-state.dtau*state.dkappa,rhs,factorization);
  

            state.dcorr_y          =d.dy;
            state.dcorr_xf         =d.dxf;
            state.dcorr_xc         =d.dxc;
            state.dcorr_tau        =d.dtau;
            state.dcorr_s          =d.ds;
            state.dcorr_kappa      =d.dkappa;
            if(isempty(state.dcorr_xf))
                state.dcorr_xf = [];
            end
                
        end
    end %End of the selection of linear solver

    %this resolves the matlab quirk that does not allow adding 
    %[] to an empty matrix
    if(isempty(state.dxf))
        state.dxf = [];
    end
   
    %Count the solution
    state.kkt_solves = state.kkt_solves + 1; 

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
        ya         = state.y        + state.a_affine*state.dy; 
        xaf        = state.xf       + state.a_affine*state.dxf;
        xca        = state.xc       + state.a_affine*state.dxc;
        taua       = state.tau      + state.a_affine*state.dtau;
        sa         = state.s        + state.a_affine*state.ds;
        kappaa     = state.kappa    + state.a_affine*state.dkappa;


        if(pars.solve_second_order)
           ya         = ya         + 0.5*state.a_affine^2*state.dcorr_y; 
           xaf        = xaf        + 0.5*state.a_affine^2*state.dcorr_xf;
           xca        = xca        + 0.5*state.a_affine^2*state.dcorr_xc;
           taua       = taua       + 0.5*state.a_affine^2*state.dcorr_tau;
           sa         = sa         + 0.5*state.a_affine^2*state.dcorr_s;
           kappaa     = kappaa     + 0.5*state.a_affine^2*state.dcorr_kappa;       
        end

        %Check if the present point is primal dual feasible
        p_feas     = eval_primal_feas(problem,xca);
        if(p_feas)
            d_feas     = eval_dual_feas(problem,sa);
            if(d_feas)
        
                %Calculate mu and the duality gap at the trial point
                dga        = xca'*sa+kappaa*taua;
                mua        = dga/(state.nu+1);
               
                %Compute s+mug(x)
                state.g     = eval_grad(problem,xca); %we can save this with a product..
                state.temp1 = +mua*state.g;
                state.temp1 = sa+state.temp1;

                %Choose the centrality measure to use
                if(pars.centrality_measure==inf)  %||s+mu g(x)||_inf/||g||_inf
                    state.cent = norm(state.temp1,inf)/norm(state.g,inf);
                    state.cent = max(state.cent,abs(taua*kappaa-mua));
                else
                    %Evaluate the centrality measure ||sa+mua*g(xa)||_H^{-1}(xa)
                    H           = eval_hessian(problem,xca);
                   
                    state.temp2 = H\state.temp1;
                    state.cent  = sqrt(state.temp1'*state.temp2);
                end
                    

                if(pars.print>3) fprintf('Bk %i ||s+mg||_H^- at affine backtrack %g,  mua*theta %g \n',b_iter, state.cent,mua*pars.theta); end

                if(state.cent <= mua*pars.theta) %Stop when this is satisfied
                    break;
                end
            else
                %not dual infeasible 
                if(pars.print>3) fprintf('Bk %i Dual infeasible at affine backtrack \n',b_iter); end
            end
        else
            %not primal infeasible 
                if(pars.print>3) fprintf('Bk %i Primal infeasible at affine backtrack \n',b_iter); end
        end 
        state.a_affine  = state.a_affine*pars.backtrack_affine_constant;
    end %End of backtrack loop

    %If the maximum number of iterates was reached report an error and exit
    if(~p_feas || ~d_feas || state.cent > mua*pars.theta)
        fprintf('Backtracking line search failed is backtrack_affine_constant too large?\n');
        state.exit_reason = 'affine backtrack line search fail';
        break;
    end

    %If the line-search succeeded take the step
    state.y     = ya;
    state.xf    = xaf;
    state.xc    = xca;
    state.tau   = taua;
    state.s     = sa;
    state.kappa = kappaa;
    state.mu    = mua;

    clear 'ya' 'xaf' 'xac' 'taua' 'sa' 'kappaa' 'mua'

    %Calculate the residuals
    state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
    state.d_res         = -problem.c*state.tau+problem.A'*state.y;
    state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
    if(~isempty(state.xf)) %If state.xf is empty then the c'x would result in an empty matrix
        cfxf = problem.c(1:nf)'*state.xf;
    else
        cfxf = 0;
    end  
    state.g_res         = - problem.b'*state.y  +cfxf+problem.c(nf+1:n)'*state.xc + state.kappa; 
    

    %          Centering process
    %----------------------------------------------------------
    if(pars.print>1) fprintf('||s+mg||_H^- before centering : %g\n',state.cent); end

    for c_iter = 1:pars.max_c_iter
        state.c_iter = c_iter;
        state.centering_iterations = state.centering_iterations +1;
        %Evaluate the gradient and hessian
        state.g = eval_grad(problem,state.xc);
        H       = eval_hessian(problem,state.xc);

        %Evaluate the centering condition
        state.temp1 = state.s+state.mu*state.g;
        state.temp2 = H\state.temp1;
        state.cent  = sqrt(state.temp1'*state.temp2);
 
        %Shorthand 
        m = problem.m;
        n = problem.n;
        nc = problem.n_constrained;
        nf = problem.n_free;


    %Solve the linear system 
    if(strcmp(pars.linear_solver,'umfpack')==1)
        %Assemble the matrix
        Sc  = [sparse(nf,nc);-speye(nc)];
        K5= [[sparse(m,m) , problem.A                           ,-problem.b     ,sparse(m,nc)  ,sparse(m,1)];...
             [-problem.A' , sparse(n,n)                         ,problem.c      ,Sc            ,sparse(n,1)];...
             [problem.b'  , -problem.c'                         ,sparse(1,1)    ,sparse(1,nc)  ,-1         ];...
             [sparse(nc,m), [sparse(nc,nf) ,sparse(state.mu*H)] ,sparse(nc,1)   ,speye(nc,nc)  ,sparse(nc,1)];...
             [sparse(1,m) , sparse(1,n)                         ,state.kappa    ,sparse(1,nc)  ,state.tau   ]];
         
              
        %Build the rhs
        rhs                 = zeros(n+m+nc+2,1);
        %XXX A little wasteful temporary solution
        r5                  = -state.s-state.mu*state.g;
        rhs(m+n+2:m+n+1+nc) = r5;
        rhs(m+n+nc+2)       = state.mu - state.tau*state.kappa;

        d       = K5\rhs; 
        state.dy       = d(1:m);
        state.dxf      = d(m+1:m+1+nf);
        state.dxc      = d(m+nf+1:m+n);
        state.dtau     = d(m+n+1);
        state.ds       = d(m+n+2:m+n+1+nc);
        state.dkappa   = d(m+n+nc+2);

    elseif(strcmp(pars.linear_solver,'mixed')==1)

        %Build the rhs
        r1             = state.p_res-(problem.b*state.tau-problem.A*[state.xf;state.xc]);
        r2             = state.d_res-(-problem.c*state.tau+problem.A'*state.y);
        r2(nf+1:n)     = r2(nf+1:n) - state.s;
        if(~isempty(state.xf)) %If state.xf is empty then the c'x would result in an empty matrix
            cfxf = problem.c(1:nf)'*state.xf;
        else
            cfxf = 0;
        end  
        r3             = state.g_res -(- problem.b'*state.y  +cfxf+problem.c(nf+1:n)'*state.xc + state.kappa);
       
        %From coneopt this is backwards
        r4             = state.mu - state.tau*state.kappa;
        r5             = -state.s-state.mu*state.g;

        d              = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,r1,r2,r3,r4,r5,[]);
        state.dy       = d.dy;
        state.dxf      = d.dxf;
        state.dxc      = d.dxc;
        state.dtau     = d.dtau;
        state.ds       = d.ds;
        state.dkappa   = d.dkappa;
       
        %If we are debugging and want to see calculate the residuals and print
        if(pars.print > 2)
            n_res_1 = norm(problem.A*[state.dxf;state.dxc]-state.dtau*problem.b-r1);
            n_res_2 = norm(-problem.A'*state.dy + state.dtau*problem.c - [zeros(nf,1);state.ds] - r2);
            n_res_3 = norm(problem.b'*state.dy-problem.c'*state.dxc -state.dkappa-r3);
            n_res_5 = norm(state.mu*H*state.dxc+state.ds-r5);
            n_res_4 = norm(state.kappa*state.dtau+state.tau*state.dkappa-r4);
            fprintf('Residuals of mixed C solve r1 %g, r2 %g, r3 %g, r5 %g, r4 %g \n',n_res_1,n_res_2,n_res_3,n_res_5,n_res_4); 
        end
    end %End of linear solver selection 

    %this resolves the matlab quirk that does not allow adding 
    %[] to an empty matrix
    if(isempty(state.dxf))
        state.dxf = [];
    end
            
    %Count the solution
    state.kkt_solves = state.kkt_solves + 1; 

    %Calculate ||dx||_H(X)
    %state.n_x_H = sqrt(state.dxc'*rhs(m+n+2:m+n+1+nc)/state.mu); %This is a cheaper way to do it rhs(m+n+2:m+n+1+nc) = -s-mu*g 
    state.n_x_H  = sqrt(state.dxc'*r5/state.mu);

    if(pars.print > 2) fprintf('\t ||dx||_H before cent backtrack : %g\n',state.n_x_H); end

    %Centering backtrack 
    %----------------------------------------------------------------------
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
 
    if(pars.print > 2) fprintf('\t ||s+mg||_H^- before cent backtrack : %g\n',state.cent); end
    
    %Backtracking iteration
    c_backtrack_iter = 0;
    for c_backtrack_iter = 1:pars.max_centering_backtrack_iter
        state.c_backtrack_iter = c_backtrack_iter;

        %Calculate the trial point
        xca        = state.xc        + state.a_cent*state.dxc;
        taua       = state.tau       + state.a_cent*state.dtau;
        sa         = state.s         + state.a_cent*state.ds;
        kappaa     = state.kappa     + state.a_cent*state.dkappa;
                  
        %Evaluate the centrality measure ||sa+mua*g(xa)||_H^{-1}(xa)
        H           = eval_hessian(problem,xca);
        state.g     = eval_grad(problem,xca); %we can save this with a product..
        
        %Observe that mu is not the value of xca'sa + taua*kappaa but the value 
        % of mu before centering starts
        state.temp1 = sa+state.mu*state.g;
        state.temp2 = H\state.temp1;
        cent        = sqrt(state.temp1'*state.temp2);
        
        if(pars.print > 3) fprintf('\t \t ||s+mg||_H^- at backtrack : %g, %i \n',state.cent,c_backtrack_iter); end
        
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
    %XXX
    %state.mu    = mua;
    
    if(pars.print > 2) fprintf('\t ||s+mg||_H^- at end of backtrack %g backtracks %i centering iter %i\n',state.cent,c_backtrack_iter,c_iter); end
    
    if(state.n_x_H <= pars.beta) %Stop centering when the centering measure ||dx||_H is smaller than beta
        break;
    end

  end %End of centering iteration
  if(pars.print > 1) fprintf('||s+mg||_H^- at end of centering : %g, main iter %i\n',state.cent,m_iter); end
 
  %Check if the centering iteration failed   
  if(state.n_x_H > pars.beta)
        fprintf('Centering failed \n');
        state.exit_reason = 'Centering fail';
        break; 
  end

    %Calculate the residuals
    state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
    state.d_res         = -problem.c*state.tau+problem.A'*state.y;
    state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
    if(~isempty(state.xf)) %If state.xf is empty then the c'x would result in an empty matrix
        cfxf = problem.c(1:nf)'*state.xf;
    else
        cfxf = 0;
    end  
    state.g_res         = - problem.b'*state.y  +cfxf+problem.c(nf+1:n)'*state.xc + state.kappa; 
    
    %Calculate the residual norms
    state.n_p_res       = norm(state.p_res,'inf')/state.rel_p_res;
    state.n_d_res       = norm(state.d_res,'inf')/state.rel_d_res;
    state.n_g_res       = abs(state.g_res)/state.rel_g_res;

    state.ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
    state.bty           = problem.b'*state.y;
    state.relative_gap  = abs( state.ctx - state.bty )/( state.tau + abs(state.bty) );

    %Print the iteration log
    %iter centering iter, 
   if(pars.print>0) 
    fprintf('%2i  %2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e\n',...
                         state.m_iter,state.c_iter,...
                         state.a_affine,state.mu,...
                         state.tau,state.kappa,...
                         state.n_p_res,state.n_d_res,...
                         state.relative_gap,state.n_g_res);
    end
    %Evaluate the stopping criteria
    if(state.n_p_res < pars.stop_primal && state.n_d_res < pars.stop_dual)
        if(state.relative_gap<pars.stop_gap)
            state.exit_reason = 'Optimal';
            break;
        elseif(state.n_g_res < pars.stop_gap && state.tau<pars.stop_tau_kappa*1.e-2*max(1,state.kappa))
            %In this case it is infeasible, try to detect if it is primal or dual infeasible
            if(state.ctx < -eps && state.bty < -eps)
                state.exit_reason  = 'Dual Infeasible';
            elseif(state.ctx > eps && state.bty > eps)
                state.exit_reason  = 'Primal Infeasible';
            end
            state.exit_reason      = 'Infeasible but undescernible';
            break;
        end     
    end
   if(state.mu < state.mu0*pars.stop_mu*1.e-2&&state.tau<1.e-2*min(1,state.kappa))
        state.exit_reason = 'Ill Posed';
        break;
   end
end %End of main loop
 
    %Print the final message
if(pars.print>0) 
    fprintf('==========================================================================================\n');
    fprintf('Exit because %s \n Iterations %i\n Centering Iterations %i\n Number of KKT Solves %i\n', state.exit_reason, state.m_iter,...
            state.centering_iterations,state.kkt_solves);
end


