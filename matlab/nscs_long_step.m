%NSCS Long step matlab minimal version
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
    pars.backtrack_affine_constant = 0.8;   %Affine backtracking constant

    %XXX: changed from 0.98 for gp testing
    pars.eta        = 0.5;                %Multiple of step to the boundary
    pars.use_nesterov_todd_centering = false; %Use centering points for symmetric cones
    pars.stop_primal= 1e-5;                 %Stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual  = 1e-5;
    pars.stop_gap   = 1e-5;
    pars.stop_mu    = 1e-5;
    pars.stop_tau_kappa = 1.e-7;
    pars.solve_second_order = true;

    pars.print      = 1;                     %Level of verbosity from 0 to 11
    %Regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 20;

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
    fprintf('Experimental NSCS long step code\n');
    fprintf(' Problem size (%i,%i) nnz(A): %i \n',problem.m,problem.n,nnz(problem.A));
    fprintf(' Free: %i, Positive %i, SOCP cones %i, Matrix %i, Exponential Cones %i\n',...
        problem.n_free,problem.n_pos,problem.n_soc_cones,...
        problem.n_sdp_cones,problem.n_exp_cones);
    fprintf('==========================================================================\n');
    fprintf('%2s  %6s   %6s   %6s     %6s     %6s       %6s       %6s     %6s     %6s\n',...
                         'it','sigma',...
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

%Sanity check verify that the dual slack is feasible
if(~eval_dual_feas(problem,state.s))
    fprintf('Error, initial dual slack not feasible');
    return;
end


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
    fprintf('%2i  %3.3e   %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e\n',...
                     state.m_iter,0,...
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
    state.g        = eval_grad(problem,state.xc);
   
    %Shorthand 
    m = problem.m;
    n = problem.n;
    nc = problem.n_constrained;
    nf = problem.n_free;
    
    %Calculate the approximate affine direction 
    %------------------------------------------------------------------------------------------
    %Call the linear solver
    [d, factorization]  = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
                         state.p_res,state.d_res,state.g_res,-state.tau*state.kappa,-state.s,[]);
 
    state.dy       = d.dy;
    state.dxf      = d.dxf;
    state.dxc      = d.dxc;
    state.dtau     = d.dtau;
    state.ds       = d.ds;
    state.dkappa   = d.dkappa;

    %this resolves the matlab quirk that does not allow adding 
    %[] to an empty matrix
    if(isempty(state.dxf))
        state.dxf = [];
    end
    
    %Count the factorization 
    state.kkt_solves = state.kkt_solves + 1; 

%    %Evaluate the second order correction
%    %And eventually the Mehrota predictor?
%    if(pars.solve_second_order)
%         
%         state.g        = eval_grad(problem,state.xc);
%        %d              = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
%        %                 zeros(size(state.p_res)),zeros(size(state.d_res)),0,-state.dtau*state.dkappa,2*state.mu*(state.xc.^(-3)).*(state.dxc.^2)+2*(-state.s-state.ds),factorization);
%        rhs             = eval_tensor(problem,state);
%        d               = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,...
%                          zeros(size(state.p_res)),zeros(size(state.d_res)),0,-state.dtau*state.dkappa,rhs,factorization);
%  
%
%        state.dcorr_y          =d.dy;
%        state.dcorr_xf         =d.dxf;
%        state.dcorr_xc         =d.dxc;
%        state.dcorr_tau        =d.dtau;
%        state.dcorr_s          =d.ds;
%        state.dcorr_kappa      =d.dkappa;
%        if(isempty(state.dcorr_xf))
%            state.dcorr_xf = [];
%        end
%            
%    end
%
    if(pars.print>2) fprintf('||dy|| %g ||dxf|| %g ||dxc|| %g ||dtau|| %g ||ds|| %g ||dkappa|| %g \n ',norm(d.dy),norm(d.dxf),norm(d.dxc),norm(d.dtau),norm(d.ds),norm(d.ds)); end
 
    clear 'K5' 'rhs' 'state.d'
    
    %Approximate tangent direction backtrack   
    %----------------------------------------------------------
    state.a_affine  = 1.0;
    %Find the largest step to the boundary of tau, kappa 
    % backtrack until the x and s points are feasible

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
        xca        = state.xc       + state.a_affine*state.dxc;
        taua       = state.tau      + state.a_affine*state.dtau;
        sa         = state.s        + state.a_affine*state.ds;
        kappaa     = state.kappa    + state.a_affine*state.dkappa;

        %Check if the present point is primal dual feasible
        p_feas     = eval_primal_feas(problem,xca);
        if(p_feas)
            d_feas     = eval_dual_feas(problem,sa);
            if(d_feas) %If primal and dual feasible
                
                %Calculate mu and the duality gap at the affine step point
                dga        = xca'*sa+kappaa*taua;
                mua        = dga/(state.nu+1);
                break;
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
    if(~p_feas || ~d_feas )
        fprintf('Backtracking line search failed is backtrack_affine_constant too large?\n');
        state.exit_reason = 'affine backtrack line search fail';
        break;
    end
    
    clear 'xac' 'taua' 'sa' 'kappaa'

    %Solve the mixed centering correcting direction
    %------------------------------------------------------------
    
    %evaluate the centering parameter using mehrotra's heuristic
    sigma = (1-state.a_affine)^3;
    
    %Calculate the correction term 
    correction_term = eval_tensor(problem,state);

    %Build the rhs
    r1             = (1-sigma)*state.p_res;
    r2             = (1-sigma)*state.d_res;
    r3             = (1-sigma)*state.g_res; 
    %Legacy from coneopt these are backwards
    r4             = sigma*state.mu - state.tau*state.kappa - (1-sigma)*state.dtau*state.dkappa;
    r5             = -state.s-sigma*state.mu*state.g;
    if pars.solve_second_order 
        r5             = r5 + (1-sigma)^3*0.5*correction_term;
    end

    
    %Call the solver be sure to re-use the factorization
    d              = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,r1,r2,r3,r4,r5,factorization);
    state.dy       = d.dy;
    state.dxf      = d.dxf;
    state.dxc      = d.dxc;
    state.dtau     = d.dtau;
    state.ds       = d.ds;
    state.dkappa   = d.dkappa; 
    %this resolves the matlab quirk that does not allow adding 
    %[] to an empty matrix
    if(isempty(state.dxf))
        state.dxf = [];
    end
   
    %If we are debugging and want to see calculate the residuals and print
    if(pars.print > 2)
        n_res_1 = norm(problem.A*[state.dxf;state.dxc]-state.dtau*problem.b-r1);
        n_res_2 = norm(-problem.A'*state.dy + state.dtau*problem.c - [zeros(nf,1);state.ds] - r2);
        n_res_3 = norm(problem.b'*state.dy-problem.c'*state.dxc -state.dkappa-r3);
        n_res_5 = norm(state.mu*H*state.dxc+state.ds-r5);
        n_res_4 = norm(state.kappa*state.dtau+state.tau*state.dkappa-r4);
        fprintf('Residuals of mixed C solve r1 %g, r2 %g, r3 %g, r5 %g, r4 %g \n',n_res_1,n_res_2,n_res_3,n_res_5,n_res_4); 
    end
   
    %Backtrack to find a feasible point 
    %----------------------------------
    state.a_affine  = 1.0;
    %Find the largest step to the boundary of tau, kappa 
    % backtrack until xc and s are feasible
    
    if(state.dtau<0)
        state.a_affine = min(state.a_affine,-state.tau/state.dtau);
    end
    if(state.dkappa<0)
        state.a_affine = min(state.a_affine,-state.kappa/state.dkappa);
    end
     
    % Backtrack loop
    b_iter = 0;
    for b_iter = 1:pars.max_affine_backtrack_iter
        state.b_iter = b_iter;

        %Evaluate the trial point
        xca        = state.xc       + state.a_affine*state.dxc;
        sa         = state.s        + state.a_affine*state.ds;

        %Check if the present point is primal dual feasible
        p_feas     = eval_primal_feas(problem,xca);
        if(p_feas)
            d_feas     = eval_dual_feas(problem,sa);
            if(d_feas) %If primal and dual feasible                
               break;
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
    if(~p_feas || ~d_feas )
        fprintf('Backtracking line search failed is backtrack_affine_constant too large?\n');
        state.exit_reason = 'affine backtrack line search fail';
        break;
    end
    
    %Take a multiple of the feasible step length just to be sure the 
    %next iterate is not too close from the boundary
    state.a_affine = state.a_affine*pars.eta;

    %Take the step
    state.y     = state.y         + state.a_affine*state.dy;
    state.xf    = state.xf        + state.a_affine*state.dxf;
    state.xc    = state.xc        + state.a_affine*state.dxc; 
    state.tau   = state.tau       + state.a_affine*state.dtau;
    state.s     = state.s         + state.a_affine*state.ds;
    state.kappa = state.kappa     + state.a_affine*state.dkappa;
    
    %Calculate mu and the duality gap at the new point
    dga              = state.xc'*state.s+state.kappa*state.tau;
    state.mu        = dga/(state.nu+1);
 
    %------------------------------------------------------------
    %Calculate the residuals
    state.p_res         =  problem.b*state.tau-problem.A*[state.xf;state.xc];
    %We just need the norm of ph_res so make it more efficient
    state.ph_res        =  problem.A*[state.xf;state.xc];
    state.d_res         = -problem.c*state.tau+problem.A'*state.y;
    state.d_res(nf+1:n) = state.d_res(nf+1:n) + state.s;
    %We just need the norm of dh_res so make it mrore efficient
    state.dh_res         = problem.A'*state.y;
    state.dh_res(nf+1:n) = state.dh_res(nf+1:n) + state.s;

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
    

    state.n_ph_res      = norm(state.ph_res,'inf')/state.rel_p_res;
    state.n_dh_res      = norm(state.dh_res,'inf')/state.rel_d_res;

    state.ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
    state.bty           = problem.b'*state.y;
    state.relative_gap  = abs( state.ctx - state.bty )/( state.tau + abs(state.bty) );
   
    %Print the iteration log
    %iter centering iter, 
   if(pars.print>0) 
    fprintf('%2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e\n',...
                         state.m_iter,sigma,...
                         state.a_affine,state.mu,...
                         state.tau,state.kappa,...
                         state.n_p_res,state.n_d_res,...
                         state.relative_gap,state.n_g_res);
    end
 
    %Evaluate the stopping criteria
    %---------------------------------------------------------
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


