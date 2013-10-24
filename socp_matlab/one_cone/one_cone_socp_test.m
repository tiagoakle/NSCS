clear all
close all

%This is a Mehrotra predictor-corrector implementation for a single cone

%SOCP solver for one cone
n = 10; %Size of the cone
m = 5; %Number of costraints
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
A = randn(m,n);
b = A*x_feas;
s = rand(n,1);
c = A'*ones(m,1) + s_feas;
fprintf('Generated a random feasible SOCP with %i constraints and %i variables\n',m,n);
clear x_feas s_feas

%-------------------
%Choose an initial x and s
x = randn(n,1);
x(1) = norm(x(2:end))+rand(1);

%Generate a feasible dual point
s    = randn(n,1);
s(1) = norm(s(2:end))+rand(1);

tau = 1;
kappa = 1;
fprintf('Selected random x,s \n');
y  = ones(m,1);

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

fprintf('%2i a %3.3g pr %3.6g dr %3.6g gr %3.3g mu %3.6g \n',0,0,nrp,nrd,nrg,mu);

%Algorithm parameters
max_iter = 40;

%--------------------------------------
%Main SOCP solver iteration
for iter=1:max_iter
   
    %Evaluate the hessian and scaled variables
    [g,H,W,Wi,lambda] = eval_scaled_variables(n,x,s);
     
    %Build the matrix
    K5= [[sparse(m,m) , A            , -b                 ,sparse(m,n)  ,sparse(m,1)];...
         [-A'         , sparse(n,n)  , c                  ,-speye(n)    ,sparse(n,1)];...
         [b'          , -c'          , sparse(1,1)        ,sparse(1,n)  ,-1         ];...
         [sparse(n,m) , sparse(H)    , sparse(n,1)        ,speye(n,n)   ,sparse(n,1)];...
         [sparse(1,m) , sparse(1,n)  ,    kappa/tau       ,sparse(1,n)  ,1          ]];

    %build the rhs for the affine direction
    rhs = [rp;...
           rd;... 
           rg;...
           -s;...
           -kappa];

    %Solve the system
    d = K5\rhs;
    dy     = d(1:m);
    dx     = d(m+1:n+m);
    dtau   = d(m+n+1);
    ds     = d(m+n+2:m+2*n+1);
    dkappa = d(m+2*n+2);
    
    %Calculate the maximum step length 
    a_max = max_step_size(tau,kappa,dx,ds,dtau,dkappa,W,Wi,lambda);
    a_max = min(1,a_max);

    mu_aff = (x+a_max*dx)'*(s+a_max*ds)+(tau+a_max*dtau)*(kappa+a_max*dkappa);
    mu_aff = mu_aff/(nu+1);
    
    %Now calculate the centering direction
    %-------------------------------------
    sigma = (1-a_max)^3;
   
   %Build the second order correction term, this should have a better form
   ws       = Wi*ds;
   wx       = W*dx;
   wsdwx    = [ws'*wx;ws(1)*wx(2:n)+wx(1)*ws(2:n)];
   li       = 1/(lambda(1)^2-norm(lambda(2:n))^2)*[lambda(1);-lambda(2:n)];
   lidsdx   = [li'*wsdwx;li(1)*wsdwx(2:n)+wsdwx(1)*li(2:n)];
   so       = W*lidsdx;

%    % %build the rhs for the centering direction
%    rhs = [(1-sigma)*rp;... 
%           (1-sigma)*rd;...
%           (1-sigma)*rg;...
%           -s-sigma*mu*g-so;...
%           -kappa-sigma*mu/tau-dtau*dkappa/tau];
    %XXXCheck this  

    % %build the rhs for the centering direction with no second order correction
    rhs = [(1-sigma)*rp;... 
           (1-sigma)*rd;...
           (1-sigma)*rg;...
           %-s-sigma*mu*g;...
           %-kappa+sigma*mu/tau];
           -s-sigma*mu*g-so;...
           -kappa+sigma*mu/tau-dtau*dkappa/tau];
       
      
    %Solve the centering direction
    d = K5\rhs;
    dy_c     = d(1:m);
    dx_c     = d(m+1:n+m);
    dtau_c   = d(m+n+1);
    ds_c     = d(m+n+2:m+2*n+1);
    dkappa_c = d(m+2*n+2);
    
    %Decrement 
    a_max = max_step_size(tau,kappa,dx_c,ds_c,dtau_c,dkappa_c,W,Wi,lambda);
    a       = a_max*0.99; 
    
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

    fprintf('%2i a %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e \n',iter,a,nrp,nrd,nrg,mu);
    if(nrp < 1.e-8 && nrd < 1.e-8 && mu <mu0*1.e-8)
        break;
    end
end %End centering


