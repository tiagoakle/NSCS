%Multi cone SOCP solver test
clear all
close all

%This is a Mehrotra predictor-corrector implementation 

%SOCP solver for multiple cones
n = 0; %Total number of variables
n_soc_cones = 2; %Total number of cones
soc_cones   = [10,10]; %Dimensions of the cones
m = 50; %Number of costraints
nu = n_soc_cones; %Complexity of the barrier
%Populate the number of variables

function [x,y,s,info] = socp_solver(A,b,c,socp_cones)
    
    
    %Deduce the problem sizes from the data 
    [m,n] = size(A);
    n_socp_ones = length(socp_cones);
    x = zeros(n,1);
    s = zeros(n,1);
    y = zeros(m,1);
    info = struct;

    %Compute the complexity parameter 
    nu = n_socp_cones;

    %Check the sizes of the variables
    if(sum(socp_cones)~=n)
       fprintf('the dimensions of A and the sizes of the cones do not match');
       return;
    end
    if(m~=size(b,1))
       fprintf('the dimensions of A and b are incompatible');
       return;
    end
   if(n~=size(c,1))
       fprintf('the dimensions of A and c are incompatible');
       return;
    end

    %-----------------------------------
    %Initialization strategy from CVXOPT
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[zeros(n,1);b];
    x   = sol(1:n);
    
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[c;zeros(m,1)];
    y   = sol(n+1:n+m);
    s   = sol(1:n);
    tau   = 1;
    kappa =1;
    clear sol
   
    %Find the indices of all the first variables 
    %for all the cones
    f_vars = [1;cumsum(n_socp_cones(1:end-1)+1];

    %Calculate the minimum shift into feasibility
    alpha = minimum_shift_into_feasilbilty(x,soc_cones);
    %Shift all primal cones into feasibility;
    x(f_vars) = x(f_vars) + alpha+1;

    %Repeat for the dual cones
    alpha = minimum_shift_into_feasilbilty(s,soc_cones);
    s(f_vars) = s(f_vars) + alpha+1;

     
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
    
    fprintf('%2i a %3.3g pr %3.6g dr %3.6g gr %3.3g mu %3.6g \n',0,0,nrp,nrd,nrg,mu);
    
    %Algorithm parameters
    max_iter = 40;
    
    %--------------------------------------
    %Main SOCP solver iteration
    for iter=1:max_iter
       
        %Evaluate the hessian and scaled variables
        [g,H,W,Wi,lambda] = eval_scaled_variables(n,x,s);
        %Build the rhs term 

        %Solve the affine scaling direction
        [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-s,-kappa,[]);
    
        %Calculate the maximum step length 
        a_max = max_step_size(tau,kappa,dx,ds,dtau,dkappa,W,Wi,lambda);
        a_max = min(1,a_max);
    
        %Now calculate the centering direction
        %-------------------------------------
        sigma = (1-a_max)^3;
    
       %Build the second order correction term, this should have a better form
        %Invert Arr(lambda)
        detl     = (lambda(1)^2-norm(lambda(2:n))^2);
        li       = [[lambda(1),-lambda(2:n)'];[-lambda(2:n),1/lambda(1)*(detl*eye(n-1)+lambda(2:n)*lambda(2:n)')]]; 
        li       = 1/detl*li;
        ws       = Wi*ds;
        wx       = W*dx;
        wsdwx    = [ws'*wx;ws(1)*wx(2:n)+wx(1)*ws(2:n)];
        lidsdx   = li*wsdwx;
        so       = W*lidsdx;
        kt_so    = dtau*dkappa/tau;
          
        %Solve the affine scaling direction
        [dy_c,dx_c,dtau_c,ds_c,dkappa_c,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,(1-sigma)*rp,...
                                                                                                (1-sigma)*rd,...
                                                                                                (1-sigma)*rg,...
                                                                                                -s-sigma*mu*g-so,...
                                                                                                -kappa+sigma*mu/tau-kt_so,...
                                                                                                slv_aug);
    
    
        %Decrement 
        a_max = max_step_size(tau,kappa,dx_c,ds_c,dtau_c,dkappa_c,W,Wi,lambda);
        a       = a_max*0.98; 
    %    a       = min(a,1); 
    
        y         = y+a*dy_c;
        x         = x+a*dx_c;
        s         = s+a*ds_c;
        tau       = tau+a*dtau_c;
        kappa     = kappa+a*dkappa_c;
    
        mu        = x'*s+tau*kappa;
        mu        = mu/(nu+1);
        
        %Feasiblity check
        x_slack = x(1)^2-norm(x(2:n))^2;
        s_slack = s(1)^2-norm(s(2:n))^2;
        if(x_slack<0||s_slack<0||x(1)<0||s(1)<0||tau<0||kappa<0)
            fprintf('Infeasible centering step \n');
            return;
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
        
        fprintf('%2i a %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',iter,a,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,residual_norm_c);
        if(nrp < 1.e-8 && nrd < 1.e-8 && mu <mu0*1.e-8)
            break;
        end
    end %End centering
end 

%Calculate the minimum alpha such that x+alpha e \in K for all primal and dual cones
%returns a negative alpha if the point is feasible
function [alpha] = minimum_shift_into_feasilbilty(x,soc_cones)
    alpha = -1;
    n_soc_cones = length(soc_cones);
    ie = 1;
    for(k=1:n_soc_cones)
        slack = x(ie)-norm(x(ie+1:ie+soc_cones(k)-1));
        alpha = max([-slack,alpha]);
        ie = ie+soc_cones(k);
    end
end

%Evaluates the scaled variables u and lambda 
%u is the square root of the scaling point w st H(w)x = s.
function [u,l,uju] = eval_scaled_variables(x,s,socp_cones)
end
