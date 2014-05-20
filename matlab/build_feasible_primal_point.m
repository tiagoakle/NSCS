function x0c = build_feasible_primal_point(problem);
        x0c = zeros(problem.n_free+problem.n_pos+3*problem.n_exp_cones,1);
        %Leave n_free zeros
        %set n_pos to one
        x0c(problem.n_free+1:problem.n_pos+problem.n_free) = ones(problem.n_pos,1);
        %Set 3*n_exp_cones to the cone center
        %This point satisfies x = -g(x)
        initial_ex=[-0.827838387734678;...
                     1.290927717327200;...
                     0.805102001257750];
        x0c(problem.n_pos+problem.n_free+1:problem.n_free+problem.n_pos+problem.n_exp_cones) = ones(problem.n_exp_cones,1)*initial_ex(1);
        x0c(problem.n_pos+problem.n_free+problem.n_exp_cones+1:problem.n_free+problem.n_pos+problem.n_exp_cones*2) = ones(problem.n_exp_cones,1)*initial_ex(2);
        x0c(problem.n_pos+problem.n_free+2*problem.n_exp_cones+1:problem.n_free+problem.n_pos+problem.n_exp_cones*3) = ones(problem.n_exp_cones,1)*initial_ex(3);
end

