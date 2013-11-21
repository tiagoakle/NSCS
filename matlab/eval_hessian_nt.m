function [H] = eval_hessian_nt(problem,xc,s,w,mu)
    %Evaluates the block hessian 
    %[1\muH_l 0 ]
    %[0      H_e]
    % where H_l is the hessian at the nt scaling point of the LP cones
    % and H_e is the hessian of the exponential cone part

    %This is wrong for socp and sdp non zero
    nnz = problem.n_pos+sum(problem.soc_cones.^2)+sum(problem.sdp_cones.^4)+...
          problem.n_exp_cones*9+problem.n_power_cones*9;

    pos_indices_1 = [];
    pos_values  = [];
    exp_indices_1 = [];
    exp_indices_2 = [];
    exp_values  = [];

    if(problem.n_pos>0)
        pos_indices_1 = [1:problem.n_pos]';
        pos_values   = 1/mu*s(1:problem.n_pos)./(xc(1:problem.n_pos));
    end
    
    %Build the block diagonal matrix in the sparity inducing form
    %Each block H(w) is of the form D + uu'. Build instead the matrix
    % [-1 u']
    % [ u D ]
    %which has exactly p+1+2p non zeros rather than p^2.

    %The hessian of the barrier is of the form 
    %d = 1/(u'Ju)
    %2d^2Juu'J -dJ
    %The diagonal is dJ
    %The outer product u = sqrt(2)d

    %The system 
    %[-1 u'][g] = [0]
    %[u  D ][x]   [b] 
    %Has the same solution x as the original Hx=b

    ie = problem.n_pos+1; %First socp variable
    ib = problem.n_pos+1; %First diagonal block variable ix
    il = 1;               %Linear index in the coo representation
    soc_indicesi = zeros(3*sum(problem.soc_cones))+problem.n_soc_cones;
    soc_indicesj = zeros(3*sum(problem.soc_cones))+problem.n_soc_cones;
    soc_vals     = zeros(3*sum(problem.soc_cones))+problem.n_soc_cones;
    if(problem.n_soc_cones>0)
       for(k=1:problem.n_soc_cones)
            %Calculate the determinant uJu
            d = xc(ie)^2-norm(xc(ie+1:ie+problem.soc_cones(k)-1))^2;
            d = 1/d;
            %Append outer product
            %initial -1
            soc_vals(il)     = -1;
            soc_indicesi(il) = ib;
            soc_indicesj(il) = ib;
            il = il+1;
            %Vertical element 
                %First vertical
                soc_vals(il)                             =  sqrt(2)*d*w(ie);
                %Remaining vertical 
                soc_vals(il+1:il+problem.soc_cones(k)-1) = -sqrt(2)*d*w(ie+1:ie+problem.soc_cones(k)-1);
                soc_indicesi(il:il+problem.soc_cones(k)-1) = [ib+1:ib+problem.soc_cones(k)]';
                soc_indicesj(il:il+problem.soc_cones(k)-1) = ib;
                il = il+problem.soc_cones(k);
            %Horizontal element
                %First horizontal
                soc_vals(il)                             =  sqrt(2)*d*w(ie);
                %Remaining horizontal
                soc_vals(il+1:il+problem.soc_cones(k)-1) = -sqrt(2)*d*w(ie+1:ie+problem.soc_cones(k)-1);
                soc_indicesi(il:il+problem.soc_cones(k)-1) = ib;
                soc_indicesj(il:il+problem.soc_cones(k)-1) = [ib+1:ib+problem.soc_cones(k)]';
                il = il+problem.soc_cones(k);  
            %Diagonal element
                %First diagonal element
                soc_vals(il) = -d;
                soc_indicesi(il) = ib+1;
                soc_indicesj(il) = ib+1;
                il = il+1;
                %p-1 diagonal elements
                soc_vals(il:il+problem.soc_cones(k)-2)     = -d;
                soc_indicesi(il:il+problem.soc_cones(k)-2) = [ib+2:ib+problem.soc_cones(k)]';
                soc_indicesj(il:il+problem.soc_cones(k)-2) = [ib+2:ib+problem.soc_cones(k)]';
                
            ie = ie+problem.soc_cones(k);
            ib = ib+problem.soc_cones(k)+1;
       end
    end

    if(problem.n_exp_cones>0)
        %Index of the first exponential variable
        i_e = problem.n_pos+sum(problem.soc_cones)+problem.n_soc_cones+sum(problem.sdp_cones.^2)+1; 
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
