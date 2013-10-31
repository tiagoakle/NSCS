clear all
close all

%This is a Mehrotra predictor-corrector implementation for a single cone

%SOCP solver for one cone
n = 100; %Size of the cone
m = 50; %Number of costraints
nu = 1; %Complexity of the barrier

%-------
%Generate a feasible primal point
x = randn(n,1);
x(1) = norm(x(2:end))+rand(1);
x_feas = x;

%Generate a feasible dual point
s = randn(n,1);
s(1) = norm(s(2:end))+rand(1);
s_feas = s;

%Generate a feasible and bounded primal dual problem
A = sprandn(m,n,0.5);
b = A*x_feas;
s = rand(n,1);
c = A'*ones(m,1) + s_feas;
fprintf('Generated a random feasible SOCP with %i constraints and %i variables\n',m,n);
clear x_feas s_feas

%%-------------------
%%Choose an initial x and s
%x = randn(n,1);
%x(1) = norm(x(2:end))+rand(1);
%
%%Generate a feasible dual point
%s    = randn(n,1);
%s(1) = norm(s(2:end))+rand(1);
%
%tau = 1;
%kappa = 1;
%
%%fprintf('Selected random x,s \n');
%y  = ones(m,1);
%%

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

if(x(1)<=0||x(1)<=norm(x(2:n))) %if x is not feasible shift it into feasibility
    x(1) = norm(x(2:n)) + 1;
end

if(s(1)<=0||s(1)<=norm(s(2:n))) %if s is not feasible shift it into feasibility
    s(1) = norm(s(2:n)) + 1;
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

    mu_aff = (x+a_max*dx)'*(s+a_max*ds)+(tau+a_max*dtau)*(kappa+a_max*dkappa);
    mu_aff = mu_aff/(nu+1);
    
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


