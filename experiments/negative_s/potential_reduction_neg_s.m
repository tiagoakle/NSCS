%Potential reduction algorithm in the homogeneous form
%This solver allows the s variables to become negative as is proceeds, 
% we are testing YYY hypothesis that the algorithm would still converge making 
% the dual slacks feasible in the long run.
function [x,y,s,tau,kappa,info] = potential_reduction_neg_s(A,b,c)

    warning('off','MATLAB:nearlySingularMatrix');

    fprintf('Minimal mehrotra type solver \n');
    %-----------------------------------------------------------------
    % Minimal Mehrotra LP predictor corrector
    %------------------------------------------------------------------
    [m,n] = size(A); 
    nu    = n;

    %The aggressiveness of rho matters here
    rho   = (nu+1)+nu*nu*sqrt(nu+1);

    fprintf('Problem size %i constraints %i variables \n',m,n);
   
    if(false)
      fprintf('CVX Opt initialization')
      %-----------------------------------
      %Initialization strategy from CVXOPT
      K2  = [[speye(n), A'];[A,sparse(m,m)]];
      sol = K2\[zeros(n,1);b];
      x   = sol(1:n);
      
      K2  = [[speye(n), A'];[A,sparse(m,m)]];
      sol = K2\[c;zeros(m,1)];
      y   = sol(n+1:n+m);
      s   = sol(1:n);
    else 
      fprintf('Badly centered initialization')
      x = rand(n,1);
      s = x;
      y = randn(m,1);
      x = x.^3;
      %%XXX:RANDOM INIT
      %x   = ones(n,1);
      %s   = ones(n,1);
      %y   = zeros(m,1);
    end 
   
    tau   = 1;
    kappa =1;
    clear sol
    
    if(min(x)<0) %if x is not feasible shift it into feasibility
        a = min(x);
        x = (1-a)*ones(n,1)+x;
    end
    
    if(min(s)<0) %if s is not feasible shift it into feasibility
        a = min(s);
        s = (1-a)*ones(n,1)+s;
    end
    
    %------Calculate the initial mu
    mu = (s'*x+tau*kappa)/(nu+1);
    mu0 = mu;
    
    %--------------------
    %Evaluate the residuals
    rp = tau*b-A*x;
    rd = -tau*c+A'*y+s;
    rg = kappa-b'*y+c'*x;
    
    %Evaluate the residual norms
    nrp = norm(rp);
    nrd = norm(rd);
    nrg = norm(rg);
    
    % Save the original norms
    nrp0 = nrp;
    nrd0 = nrd;
    nrg0 = nrg;
    
    gap  = c'*x-b'*y;
    fprintf('%2i a       %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',0,nan,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,nan);
    
    %Algorithm parameters
    max_iter = 1000;
    sigma    = (nu+1)/rho;
    
    %--------------------------------------
    %Main solver iteration
    for iter=1:max_iter
       
        %Evaluate the hessian and scaled variables
        %Evaluate the scaling point 
        H = mu*diag(sparse((1./x).*(1./x)));
       
          
        %Solve the affine scaling direction
        [dy_c,dx_c,dtau_c,ds_c,dkappa_c,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,...
                                                                                                (1-sigma)*rp,...
                                                                                                (1-sigma)*rd,...
                                                                                                (1-sigma)*rg,...
                                                                                                -s    +sigma*mu*1./x,...
                                                                                                -kappa+sigma*mu*1/tau,...
                                                                                                []);
        %Calculate the newton decrement
        lambda = norm([(dx_c.*1./x);sqrt(dtau_c*kappa/tau*dtau_c)]);
        
        ratios = [-x./dx_c;-tau/dtau_c;-kappa/dkappa_c;1];
        a_max  = min(ratios(find(ratios>0)));
        
        %Calculate the step size a la nesterov
        a_step = 1/(1+lambda);
        if(a_step>a_max)
            fprintf('Error, nesterov step larger than primal feasible step\n');
            return;
        end
          
        x_old = x;
        y_old = y;
        s_old = s;
        t_old = tau;
        k_old = kappa;
        
        %Take the step
        a = a_step;
     
        y         = y+a*dy_c;
        x         = x+a*dx_c;
        s         = s+a*ds_c;
        tau       = tau+a*dtau_c;
        kappa     = kappa+a*dkappa_c;
    
        mu        = x'*s+tau*kappa;
        mu        = mu/(nu+1);
        
        %Feasiblity check
        x_slack = min(x);
        if(x_slack<0||tau<0||kappa<0)
            fprintf('Infeasible centering step \n');
            return;
        end
        s_slack = min(s);
        if(s_slack<0)
            fprintf('Duals are infeasible \n');
        end 
    
        %--------------------
        %Evaluate the residuals
        rp = tau*b-A*x;
        rd = -tau*c+A'*y+s;
        rg = kappa-b'*y+c'*x;
    
        nrp   = norm(rp);
        nrd   = norm(rd);
        nrg   = norm(rg);
        
        gap   = (mu*(nu+1)-tau*kappa)/tau;
        
        %Measure the centrality
        centrality = norm([x.*s-mu;tau*kappa-mu]);
        fprintf('%2i a %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e lda %3.3e cent %3.3e mins %3.3e\n',...
                iter,a,nrp/(nrp0),nrd/(nrd0),nrg,mu,gap,kappa/tau,lambda,centrality,min(s)/tau);

        if(nrp/nrp0 < 1.e-8 && nrd/nrd0 < 1.e-8 && mu/mu0 < 1.e-8 && min(s)/tau > - 1.e-10)
            break;
        end
    end %end main loop

    x = x/tau;
    s = s/tau;
    y = y/tau;
    info = struct;
    info.iter = iter;
end