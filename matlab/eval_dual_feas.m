function feas = eval_dual_feas(problem,xc)
    feas = 1;

    if(problem.n_pos>0)
        min_pos_ort = min(xc(1:problem.n_pos));
        if(min_pos_ort<=0) 
            feas = 0;
            return;
        end
    end

    if(problem.n_soc_cones>0)
        i_e=problem.n_pos+1;
        for(k=1:problem.n_soc_cones)
            if(xc(i_e)<norm(xc(i_e+1:i_e+problem.soc_cones(k)-1)))
                feas = 0;
                return;
            end
            i_e = i_e+problem.soc_cones(k);
        end
    end

    if(problem.n_exp_cones>0)
        %Index of the first exponential variable
        i_e = problem.n_pos+sum(problem.soc_cones)+sum(problem.sdp_cones.^2)+1;
        x3    = -xc(i_e:i_e+problem.n_exp_cones-1);
        x2    = exp(1)*xc(i_e+problem.n_exp_cones:i_e+2*problem.n_exp_cones-1);
        x1    = -xc(i_e+2*problem.n_exp_cones:i_e+3*problem.n_exp_cones-1);
        logx2 = log(x2);
        logx3 = log(x3);
        tmp1  = logx2-logx3;
        psi   = x3.*tmp1 - x1;

        if(min(psi)<0)
            feas = 0;
            return;
        end
        if(min(x2)<0)
            feas = 0;
            return;
        end
        if(min(x3)<0)
            feas = 0;
            return;
        end
    end
  
 
end

