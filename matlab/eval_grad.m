function [g] = eval_grad(problem,xc)
    %Evaluates the gradient of the primal barrier at the point xc
    %the barrier is determined by the problem defined in problem.
    
    %Allocate space for the gradient
    g = zeros(problem.n_constrained,1);
    %Evaluate the positive orthant variables
    g(1:problem.n_pos) = -1./xc(1:problem.n_pos);
    %Shift the index 
    ix = problem.n_pos+1;

    %Evaluate the exponential cone barrier gradient.
    %The entries are organized as 
    %in the pattern [x1s, x2s, x3s]
    if(problem.n_exp_cones>0)
       x1       = xc(ix:ix+problem.n_exp_cones-1);
       x2       = xc(ix+problem.n_exp_cones:ix+2*problem.n_exp_cones-1);
       x3       = xc(ix+2*problem.n_exp_cones:ix+3*problem.n_exp_cones-1);

       logx2 = log(x2);
       logx3 = log(x3);
       x2m1  = 1./x2;
       x3m1  = 1./x3;
       tmp1  = logx2-logx3;
       psi   = x3.*tmp1 - x1;
       psim1 = 1./psi;
       xi    = tmp1 - 1;


       tmp2  = [ psim1;...
            -x2m1.*(x3.*psim1 + 1);...
            -xi.*psim1 - x3m1];
       g(ix:ix+3*problem.n_exp_cones-1) = tmp2(:);
    end

end
