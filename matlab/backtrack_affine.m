function [state] = backtrack_affine(state,pars,problem)
    % Backtrack loop
    b_iter = 0;
    p_feas = false;
    d_feas = false;
    
    state.failed = false;

    for b_iter = 1:pars.max_affine_backtrack_iter
        state.b_iter = b_iter;
    
        %Evaluate the trial point
        xca        = state.xc       + state.a_affine*state.dxc;
        taua       = state.tau      + state.a_affine*state.dtau;
        sa         = state.s        + state.a_affine*state.ds;
        kappaa     = state.kappa    + state.a_affine*state.dkappa;
        
        %Calculate mu and the duality gap at the affine step point
        dga        = xca'*sa+kappaa*taua;
        mua        = dga/(state.nu+1);
    
        %Check if the present point is primal dual feasible
        p_feas     = eval_primal_feas(problem,xca);
        if(p_feas)
            d_feas     = eval_dual_feas(problem,sa);
            if(d_feas) %If primal and dual feasible
                cent = eval_small_neigh(problem,xca,sa,mua);
                if(cent <pars.neigh)
                    break;
                else
                    if(pars.print>3) fprintf('Bk %i Neighborhood volation at affine\n',b_iter); end
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
    if(~p_feas || ~d_feas || cent>pars.neigh)
        if(pars.print>0)
            fprintf('Backtracking line search failed is backtrack_affine_constant too large?\n');
        end
        state.exit_reason = 'affine backtrack line search fail';
        state.fail      = true;
        return;
    end
          
end
