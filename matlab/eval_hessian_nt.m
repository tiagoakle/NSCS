function [H] = eval_hessian_nt(problem,xc,s,mu)
    %Evaluates the sparse block hessian 
    %[1\muH_l 0 ]
    %[0      H_e]
    % where H_l is the hessian at the nt scaling point of the LP cones
    % and H_e is the hessian of the exponential cones evaluated at xc

    %Number of non zeros
    nnz = problem.n_pos+...
          problem.n_exp_cones*9;

    pos_indices_1 = [];
    pos_values  = [];
    exp_indices_1 = [];
    exp_indices_2 = [];
    exp_values  = [];

    if(problem.n_pos>0)
        pos_indices_1 = [1:problem.n_pos]';
        pos_values   = 1/mu*s(1:problem.n_pos)./(xc(1:problem.n_pos));
    end
    
     if(problem.n_exp_cones>0)
        %Index of the first exponential variable
        i_e = problem.n_pos+1; 
        %Indices for the columns (and rows) where the entries go
        ix1   = i_e+[0:problem.n_exp_cones-1]';
        ix2   = i_e+problem.n_exp_cones+[0:problem.n_exp_cones-1]';
        ix3   = i_e+2*problem.n_exp_cones+[0:problem.n_exp_cones-1]';

        %Calculate the 6 entries of the hessian
        x1    = xc(i_e:i_e+problem.n_exp_cones-1);
        x2    = xc(i_e + problem.n_exp_cones:i_e+2*problem.n_exp_cones-1);
        x3    = xc(i_e + 2*problem.n_exp_cones:i_e+3*problem.n_exp_cones-1);

        logx2 = log(x2);
        logx3 = log(x3);
        x2m1  = 1./x2;
        x3m1  = 1./x3;
        tmp1  = logx2-logx3;
        psi   = x3.*tmp1 - x1;
        psim1 = 1./psi;
        xi    = tmp1 - 1;
        psi2  = psi.*psi;
        psim2 = psim1.*psim1;
        x2m2  = x2m1.*x2m1;
        x3m1  = 1./x3;
        x3m2  = x3m1.*x3m1;

        el11  = psim2;
        el21  = -x3.*x2m1.*psim2;
        el31  = -xi.*psim2;
        el22  = psim2.*x2m2.*(x3.*psi + x3.*x3 + psi2);
        el32  = psim2.*x2m1.*(x3.*xi - psi);
        el33  = psim2.*(x3m1.*psi + xi.*xi + psi2.*x3m2);

        exp_indices_1 =  [ix1;...
                          ix1;...
                          ix1;...
                          ix2;...
                          ix2;...
                          ix2;...
                          ix3;...
                          ix3;...
                          ix3];

        exp_indices_2 =  [ix1;...
                          ix2;...
                          ix3;...
                          ix1;...
                          ix2;...
                          ix3;...
                          ix1;...
                          ix2;...
                          ix3];

        exp_values  =  [el11;...
                        el21;...
                        el31;...
                        el21;...
                        el22;...
                        el32;...
                        el31;...
                        el32;...
                        el33];
    end
    H = sparse([pos_indices_1;exp_indices_1],[pos_indices_1;exp_indices_2],[pos_values;exp_values]);
end
