function r = eval_tensor(problem,state)
 %Evaluates the product
 %1/mu(\nabla^3f)[H^(-1)(ds+s),H^{-1}(ds+s)] - 2(s+ds)
 %For the first order tensor
 r = zeros(problem.n_constrained,1);
 temp = zeros(problem.n_constrained,1);
 temp = state.s+state.ds;
 
 %Evaluate the tensor of the barrier of the positive orthant
 if(problem.n_pos>0)
    r(1:problem.n_pos) = 2/state.mu*state.xc(1:problem.n_pos).*temp(1:problem.n_pos).^2;
 end
 
 if(problem.n_exp_cones>0)
    %Calculate the index of the first variable of the exponential cones
    i_e = problem.n_pos + sum(problem.soc_cones) + sum(problem.sdp_cones.^2)+1;
    %Extract all the x1s
    x = state.xc(i_e:i_e+problem.n_exp_cones-1);
    y = state.xc(i_e+problem.n_exp_cones:i_e+2*problem.n_exp_cones-1);
    z = state.xc(i_e+2*problem.n_exp_cones:i_e+3*problem.n_exp_cones-1);
    %Extract the parts of the rhs
    a1 = temp(i_e:i_e+problem.n_exp_cones-1);
    a2 = temp(i_e+problem.n_exp_cones:i_e+2*problem.n_exp_cones-1);
    a3 = temp(i_e+2*problem.n_exp_cones:i_e+3*problem.n_exp_cones-1);
    %Now evaluate \nabla^3(x)[H^{-1}(x)[a;b;c],H^{-1}(x)[a;b;c]]

    lyz= log(y)-log(z);
    t1 = (a2.*y + 2.*a1.*z - a3.*z - a1.*z.*lyz);
    d1 = (2.*z - x + z.*lyz);
    
    r1 = (z.*t1.^2)./d1.^2 - 2.*a1.^2.*(x - z.*lyz);
    
    r2 = (2.*z.*(a2.*y - a1.*z).*t1)./(y.*d1) - (z.^2.*t1.^2)./(y.*d1.^2) -...
         (2.*(a2.^2.*y.^2 - a1.^2.*z.^2 - a1.^2.*x.*z + 2.*a1.^2.*z.^2.*lyz + a1.*a3.*z.^2 + a1.*a2.*y.*z))./y; 
    
    r3 = a1.^2.*x - 26.*a1.^2.*z - 2.*a3.^2.*z - 4.*a1.^2.*z.*lyz.^2 +...
         ((a1.*x - a2.*y - 4.*a1.*z + a3.*z).*(3.*a1.*x - a2.*y - 14.*a1.*z + 3.*a3.*z))./d1 -...
         6.*a1.*a2.*y + 12.*a1.*a3.*z + 2.*a1.^2.*x.*lyz + 11.*a1.^2.*z.*lyz -...
         (z.*(lyz + 1).*t1.^2)./d1.^2 - 4.*a1.*a3.*z.*log(y./z);
     
     r(i_e:i_e+3*problem.n_exp_cones-1) = -1/state.mu*[r1;r2;r3];
    
 end
 %Add the -2muH(x') term
 r = r-2*temp;
end

