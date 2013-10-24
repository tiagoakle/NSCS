%Calculate the largest step size to the boundary
function a_max = max_step_size(tau,kappa,dx,ds,dtau,dkappa,W,Wi,lambda)
         %SOCP variables
         n     = size(lambda,1);
         J     = diag(sparse([1;-ones(n-1,1)]));
         Jl    = J*lambda;
         lJl   = lambda'*Jl;
         l_h   = lambda/sqrt(lJl);
         lJ_h  = Jl/sqrt(lJl);
         s_h   = Wi*ds;
         x_h   = W*dx;
         
         ljs   = lJ_h'*s_h;
         rho   = [ljs;s_h(2:n)-(ljs+s_h(1))/(l_h(1)+1)*l_h(2:n)];
         rho   = rho/sqrt(lJl);
         ljx   = lJ_h'*x_h;
         sig   = [ljx;x_h(2:n)-(ljx+x_h(1))/(l_h(1)+1)*l_h(2:n)];
         sig   = sig/sqrt(lJl);
       
         a_max_s = max([0,norm(rho(2:n))-rho(1),norm(sig(2:n))-sig(1)]);
         a_max_s = 1/a_max_s;

         r = [-tau/dtau,-kappa/dkappa];
         a_max_tk = min(r(find(r>0)));
         if(isempty(a_max_tk)) a_max_tk = inf; end
         a_max = min(a_max_s, a_max_tk);
end
