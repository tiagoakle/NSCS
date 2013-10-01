function [H] = eval_hessian(problem,xc)
    nnz = problem.n_pos+sum(problem.soc_cones.^2)+sum(problem.sdp_cones.^2)+...
          problem.n_exp_cones*9+problem.n_power_cones*3;
    H   = sparse(problem.n_constrained,problem.n_constrained,nnz);
    H(1:problem.n_pos,1:problem.n_pos) = diag(sparse(1./(xc(1:problem.n_constrained)).^2));
end
