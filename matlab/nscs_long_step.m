%NSCS Long step matlab minimal version

function [xc,xf,y,s,t,k,info] = nscs_long_step(problem,x0f,x0c,pars)
%Problem must contain the fields  
%   m
%   n_free
%   n_pos
%   n_soc_cones
%   soc_cones[]
%   n_sdp_cones
%   sdp_cones[]
%   n_exp_cones 
%   Starting point
%   x0f, x0c 
%   Problem definition in problem strucutre
%    A
%    b
%    c
    
    %Quiet the matlab warning about bad scaling
    warning('off','MATLAB:nearlySingularMatrix');

    %-------------------------------------------------
    % Validate the problem definition
    %-------------------------------------------------
    
    [ret,problem]=validate_problem_structure(problem);
    if(ret~=0)
        return;
    end 
   
    %If pars is not defined get the default 
    if ~exist('pars') 
        pars = set_default_pars_nscs_long_step()
    end
  
    %If the initial point is not set set it
    if(isempty(x0c))
        x0c = build_feasible_primal_point(problem);
    end
      
    
    %--------------------------------------------------
    % 'state' holds the present iterates and directions
    %--------------------------------------------------
    
    %State will be a global variable so that we dont
    %make copies
    state = struct;
    state.xc = x0c;
    state.xf = x0f; 
    
    %-----------------------------------------------
    %Initialize the variables that are a function of 
    % the initial point and parameters
    %----------------------------------------------- 
    [state,problem] = setup_state(state,problem,x0f,x0c);

    %Verify that the inital primal is feasible
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
        fprintf(' Free: %i, Positive %i Exponential Cones %i\n',...
            problem.n_free,problem.n_pos,problem.n_exp_cones);
        fprintf('==========================================================================\n');
        fprintf('%2s  %6s   %6s   %6s     %6s     %6s       %6s       %6s     %6s     %6s     %6s\n',...
                             'it','sigma',...
                             'a_a','mu',...
                             'tau','kap',...
                             'p_res','d_res','rel_gap','g_res','sigma');
        fprintf('==========================================================================\n');
    end
    
    %Print the iteration log
    if(pars.print > 0)
        fprintf('%2i  %3.3e   %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e  %3.3e\n',...
                         state.m_iter,0,...
                         0,state.mu,...
                         state.tau,state.kappa,...
                         state.n_p_res,state.n_d_res,state.relative_gap,state.n_g_res,state.sigma);
    end
    
    %-----------------------------------------------
    %Start of the major iteration
    %-----------------------------------------------
    
    %Shorthand symbols
    n           = problem.n;
    nc          = problem.n_constrained;
    nf          = problem.n_free;
    m           = problem.m;

    %Iteration counter
    m_iter = 0;
    %Message of exit reason
    state.exit_reason = 'Max Iter Reached';

    for m_iter = 1:pars.max_iter
        
        state.m_iter = m_iter;
        
        %Shorthand 
        m = problem.m;
        n = problem.n;
        nc = problem.n_constrained;
        nf = problem.n_free;
        
        %-------------------------------------------
        %Calculate the approximate affine direction 
        %-------------------------------------------
        
        %Evaluate the hessian and gradient
        H              = eval_hessian_nt(problem,state.xc,state.s,state.mu); 
        state.g        = eval_grad(problem,state.xc);

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
        
        %If we're on a more verbose mode then print the norms of the step directions
        if(pars.print>2) fprintf('||dy|| %g ||dxf|| %g ||dxc|| %g ||dtau|| %g ||ds|| %g ||dkappa|| %g \n ',...
                            norm(d.dy),norm(d.dxf),norm(d.dxc),norm(d.dtau),norm(d.ds),norm(d.ds)); end
         
        clear 'K5' 'rhs' 'state.d'
        
        %---------------------------------------
        %Backtrack into feasibility
        %---------------------------------------

        state.a_affine  = 1.0;
        %Find the largest step to the boundary of the symmetric cones
        state.a_affine = max_step_symmetric_cones(problem,state);
        
        state = backtrack_affine(state,pars,problem);
        if(state.failed)
            break;
        end

        %----------------------------------------------
        %Solve the mixed centering correcting direction
        %----------------------------------------------
        
        %evaluate the centering parameter using mehrotra's heuristic
        sigma = (1-state.a_affine)^3;
        state.sigma = sigma; 
        %Calculate the correction term 
        correction_term = eval_tensor(problem,state,pars);
    
        %Build the rhs
        r1             = (1-sigma)*state.p_res;
        r2             = (1-sigma)*state.d_res;
        r3             = (1-sigma)*state.g_res; 
        %Legacy from coneopt these are backwards
        r4             = sigma*state.mu - state.tau*state.kappa - (1-sigma)*state.dtau*state.dkappa;
        r5             = -state.s-sigma*state.mu*state.g;
        
        if(pars.solve_second_order)
            %Add the correction to the exponential cone part 
            i_e = problem.n_pos + 1;
            r5(i_e:i_e+3*problem.n_exp_cones-1) = r5(i_e:i_e+3*problem.n_exp_cones-1) + ...
                (1-sigma)*0.5*correction_term(i_e:i_e+3*problem.n_exp_cones-1);
            %Add the correction for the symmetric part
            r5(1:problem.n_pos) = r5(1:problem.n_pos) + 0.5*(1-sigma)*correction_term(1:problem.n_pos);
        end
     
        %Call the solver and e-use the factorization
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
            taua       = state.tau      + state.a_affine*state.dtau;
            kappaa     = state.kappa    + state.a_affine*state.dkappa;

            dga        = xca'*sa+kappaa*taua;
            mua        = dga/(state.nu+1);

            %Check if the present point is primal dual feasible
            p_feas     = eval_primal_feas(problem,xca);
            if(p_feas)
                d_feas     = eval_dual_feas(problem,sa);
                if(d_feas) %If primal and dual feasible                
                    if(eval_small_neigh(problem,xca,sa,mua)<pars.neigh)
                        break;
                    else
                        if(pars.print>3) fprintf('Bk %i Neighborhood volation \n',b_iter); end
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

        %Calculate the present value of the centrality
        %eval_small_neigh(problem,state.xc,state.s,state.mu);
     
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
        
        state.ctx           = problem.c(problem.n_free+1:problem.n)'*state.xc;
        state.bty           = problem.b'*state.y;
        state.relative_gap  = abs( state.ctx - state.bty )/( state.tau + abs(state.bty) );
       
        %Print the iteration log
        %iter centering iter, 
       if(pars.print>0) 
        fprintf('%2i  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e  %3.3e   %3.3e   %3.3e\n',...
                             state.m_iter,sigma,...
                             state.a_affine,state.mu,...
                             state.tau,state.kappa,...
                             state.n_p_res,state.n_d_res,...
                             state.relative_gap,state.n_g_res,state.sigma);
        end
     
        %Evaluate the stopping criteria
        %---------------------------------------------------------

        if(state.n_p_res < pars.stop_primal && state.n_d_res < pars.stop_dual)
            if(state.relative_gap<pars.stop_gap)
                state.exit_reason = 'Optimal';
                break;
            elseif(state.n_g_res < pars.stop_gap && state.tau<pars.stop_tau_kappa*max(1,state.kappa))
                %In this case it is infeasible, try to detect if it is primal or dual infeasible
                if(state.ctx < -eps && state.bty < -eps)
                    state.exit_reason  = 'Dual Infeasible';
                elseif(state.ctx > eps && state.bty > eps)
                    state.exit_reason  = 'Primal Infeasible';
                else
                    state.exit_reason      = 'Infeasible but undescernible';
                end
                break;
            end     
        end
       if(state.mu < state.mu0*pars.stop_mu&&state.tau<1.e-2*min(1,state.kappa))
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

    %Save the info
    info.kkt_solves = state.kkt_solves;
    info.exit_reason = state.exit_reason;
    xc = state.xc;
    xf = state.xf;
    y  = state.y;
    s  = state.s;
    t  = state.tau;
    k  = state.kappa;
end

