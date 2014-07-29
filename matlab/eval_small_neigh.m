%this function returns 1 if the present iterate belongs 
%to the neighborhood 1/mu||s+mu g(x)||x^* 

%Here mu is the complementarity of the exponential cones, g the gradient 
%of the exponential cones and s the duals of the exponential cones

function norm = eval_small_neigh(problem,xc,s,tau,kappa,mu)
    norm = 0; %Accumulator
    sc    = zeros(3,1);

    %Calculate mu
    i_e = problem.n_pos+1;
    %stranspose x fr the conical
    for j=1:problem.n_exp_cones
        %Extract the entries of x for the cone
        x1    = xc(i_e);
        x2    = xc(i_e+problem.n_exp_cones);
        x3    = xc(i_e+2*problem.n_exp_cones);

        %Extract the entries of s for the cone
        sc(1)    = s(i_e);
        sc(2)    = s(i_e+problem.n_exp_cones);
        sc(3)    = s(i_e+2*problem.n_exp_cones);
      
        
        l     = log(x2)-log(x3);
        r     = x3*l-x1; 
        yinv  = 1/x2;
        zinv  = 1/x3;
        rinv  = 1/r;
        rinvsqr = sqrt(1/r);
        rzsqr   = sqrt(r+x3);
        rzinvsqr = 1/rzsqr;

        %Gradient
        g  = [rinv;-yinv-x3*rinv*yinv;(1-l)*rinv-zinv];
 
        %Cholesky of H 
        C =  ...
        [[        rinv,                          0,                                   0];...
         [   -x3*rinv*yinv,  rzsqr*yinv*rinvsqr,                                      0];...
         [ -(l - 1)*rinv , -1*rinvsqr*rzinvsqr     , (r + 2*x3)^(1/2)*zinv*rzinvsqr   ]];

        %Calculate (s+mu g)H^-1(s+mu g)       
        v     = linsolve(C,sc/mu+g);
        v     = linsolve(C',v);
        norm  = norm + ((sc/mu+g)'*v);
        i_e   = i_e+1;

        %%Sanity checks
        %logx2 = log(x2);
        %logx3 = log(x3);
        %x2m1  = 1./x2;
        %x3m1  = 1./x3;
        %tmp1  = logx2-logx3;
        %psi   = x3.*tmp1 - x1;
        %psim1 = 1./psi;
        %xi    = tmp1 - 1;
        %psi2  = psi.*psi;
        %psim2 = psim1.*psim1;
        %x2m2  = x2m1.*x2m1;
        %x3m1  = 1./x3;
        %x3m2  = x3m1.*x3m1;

        %el11  = psim2;
        %el21  = -x3.*x2m1.*psim2;
        %el31  = -xi.*psim2;
        %el22  = psim2.*x2m2.*(x3.*psi + x3.*x3 + psi2);
        %el32  = psim2.*x2m1.*(x3.*xi - psi);
        %el33  = psim2.*(x3m1.*psi + xi.*xi + psi2.*x3m2);
        %
        %H = [[el11,el21,el31];[el21,el22,el32];[el31,el32,el33]];
        %max(max(abs(C*C'-H)))

    end 
  %  fprintf('Norm %g\n',norm);
end
