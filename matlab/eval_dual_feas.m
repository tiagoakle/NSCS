function feas = eval_dual_feas(problem,xc)
    feas = 1;
    min_pos_ort = min(xc(1:problem.n_constrained));
    if(min_pos_ort<=0) 
        feas = -1;
        return;
    end
end

