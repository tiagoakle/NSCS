function feas = eval_primal_feas(problem,xc)
    feas = 1;
    min_pos_ort = min(xc);
    if(min_pos_ort<=0) 
        feas = 0;
        return;
    end
end

