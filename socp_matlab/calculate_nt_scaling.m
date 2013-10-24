%Calculates the Nesterov-Todd scaling point for the socp problem 
%defines in problem.
function w = calculate_nt_scaling(problem,x,s)
    w = zeros(size(x));
    i = 1;
    for(k = 1:problem.n_soc_cones)
        xc = x(i:i+problem.soc_cones(k)-1);
        sc = s(i:i+problem.soc_cones(k)-1);
        J = diag(sparse([1;-ones(problem.soc_cones(k)-1,1)]));
        Jx = J*xc;
        xJx = xc'*Jx;
        Js = J*sc;
        sJs = sc'*Js;
        sx  = sc'*xc;
        w_gamma = sqrt((1+(sx/sqrt(sJs*xJx)))/2);
        w(i:i+problem.soc_cones(k)-1) = sqrt(sqrt((xJx/sJs)))*(1/(2*w_gamma))*(xc/sqrt(xJx)+1/sqrt(sJs)*J*sc);
        i   = i+problem.soc_cones(k);
    end

end

