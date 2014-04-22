
%Potential reduction algorithm in the homogeneous form
%This solver allows the s variables to become negative as is proceeds, 
% we are testing YYY hypothesis that the algorithm would still converge making 
% the dual slacks feasible in the long run, for the exponential cone case!

function [xc,xf,y,s,info] = potential_reduction_neg_s_exponential_cone(problem,x0f,x0c,pars)


    %Quiet the matlab warning about bad scaling
    warning('off','MATLAB:nearlySingularMatrix');
   
    state = struct;
    state.xc = x0c;
    state.xf = x0f;
       
    %Verify that the inital d is feasible
    if(~eval_primal_feas(problem,state.xc))
        fprintf('Error, initial primal slack not feasible');
        return;
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
    state.nu = problem.n_pos+problem.n_soc_cones; %Each soc cone is complexity one
    state.nu = state.nu + sum(problem.sdp_cones);
    state.nu = state.nu + problem.n_exp_cones*3;

    %Calculate RHO and sigma
    %---------------------------------------
    rho = (state.nu+1)^4 + sqrt(state.nu+1) 
    sigma = (state.nu+1)/rho;
    %--------------------------------------

    %Allocate the space for the working vectors
    state.temp1 = zeros(problem.n_constrained,1);
    state.temp2 = zeros(problem.n_constrained,1);
    
    %Calculate the starting point
    state.y  = ones(problem.m,1);
    state.s  = ones(problem.n_constrained,1);
    state.xf = x0f;
    state.xc = x0c;
    state.tau = 1;
    state.kappa = 1;

    %Calculate a feasible centered dual slack
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
       
    %Allocate space for the soc scaling points 
    state.w      = zeros(sum(problem.soc_cones));
    state.w_h    = zeros(sum(problem.soc_cones)); %used to calclate products with W and W_i
    state.l      = zeros(sum(problem.soc_cones)); %Scaled variable

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
        
        %Calculate the scaling points for the NT scaling 
        [w,w_h,l]=eval_scaling_points(problem,state.xc,state.s);
        state.w = w;
        state.w_h = w_h;
        state.l = l;
       
        %Shorthand 
        m = problem.m;
        n = problem.n;
        nc = problem.n_constrained;
        nf = problem.n_free;
        
        %-------------------------------------------
        %Calculate the approximate affine direction 
        %-------------------------------------------
        
        %Evaluate the hessian and gradient
        H              = eval_hessian(problem,state.xc); 
        state.g        = eval_grad(problem,state.xc);

        %Build the RHS
        r1             = (1-sigma)*state.p_res;
        r2             = (1-sigma)*state.d_res;
        r3             = (1-sigma)*state.g_res; 
        %Legacy from coneopt these are backwards
        r4             = sigma*state.mu - state.tau*state.kappa;
        r5             = -state.s-sigma*state.mu*state.g;

        d  = solve_linear_system(H,state.mu,state.kappa,state.tau,problem,pars,r1,r2,r3,r4,r5,[]);
 

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
   
   if(pars.print>2) fprintf('||dy|| %g ||dxf|| %g ||dxc|| %g ||dtau|| %g ||ds|| %g ||dkappa|| %g \n ',...
                            norm(d.dy),norm(d.dxf),norm(d.dxc),norm(d.dtau),norm(d.ds),norm(d.ds)); end
     
        clear 'K5' 'rhs' 'state.d'
        
        %---------------------------------------
        % Calculate lambda 
        %---------------------------------------
        hdx = state.dxc'*(H*[state.dxc]) ...
              + state.dtau*state.kappa/state.tau*state.dtau;
        lambda = norm(sqrt(hdx));
        %Calculate the step size
        state.a_affine = 1/(1+lambda);
     
       
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

        %Sanity check
        p_feas     = eval_primal_feas(problem,state.xc);
        if(~p_feas)
            fprintf('Primal infeasible \n!')
            return;
        end
 
     
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
            elseif(state.n_g_res < pars.stop_gap && state.tau<pars.stop_tau_kappa*max(1,state.kappa))
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
end

