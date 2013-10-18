clear all
close all

%--------------------
%Experiment parameters
n_c = 10;
%Defie some second order cones
problem = struct;
problem.soc_cones = [n_c,2*n_c,3*n_c];
problem.n_soc_cones = 3;

m = 10;
n = sum(soc_cones);
nu= 2*problem.n_soc_cones; %Parameter of the problem

%-------
%Generate a feasible primal point
x1 = randn(problem.soc_cones(1),1);
x2 = randn(problem.soc_cones(2),1);
x3 = randn(problem.soc_cones(3),1);
x1 = norm(x1(2:end))+rand(1);
x3 = norm(x2(2:end))+rand(1);
x3 = norm(x3(2:end))+rand(1);
x_feas = [x1;x2;x3];

%Generate a feasible dual point
s1 = randn(problem.soc_cones(1),1);
s2 = randn(problem.soc_cones(2),1);
s3 = randn(problem.soc_cones(3),1);
s1 = norm(s1(2:end))+rand(1);
s3 = norm(s2(2:end))+rand(1);
s3 = norm(s3(2:end))+rand(1);
s_feas = [s1;s2;s3];


%Generate a feasible and bounded primal dual problem
A = randn(m,n);
b = A*x_feas;
s = rand(n,1);
c = A'*ones(m,1) + s_feas;
fprintf('Generated a random feasible SOCP with %i constraints and %i variables\n',m,n);

%-------------------
%Choose an initial x and s
x1 = randn(problem.soc_cones(1),1);
x2 = randn(problem.soc_cones(2),1);
x3 = randn(problem.soc_cones(3),1);
x1 = norm(x1(2:end))+rand(1);
x3 = norm(x2(2:end))+rand(1);
x3 = norm(x3(2:end))+rand(1);
x  = [x1;x2;x3];

%Generate a feasible dual point
s1 = randn(problem.soc_cones(1),1);
s2 = randn(problem.soc_cones(2),1);
s3 = randn(problem.soc_cones(3),1);
s1 = norm(s1(2:end))+rand(1);
s3 = norm(s2(2:end))+rand(1);
s3 = norm(s3(2:end))+rand(1);
s  = [s1;s2;s3];
tau = 1;
kappa = 1;
fprintf('Selected random x,s \n');
y  = ones(m,1);

%------Select a target mu
mu = 0.05*(s'*x+tau*kappa)/(nu+1);

%--------------------
%Evaluate the residuals
rp = tau*b-A*x;
rd = -tau*c+A'*y+s;
rg = kappa-b'*y+c'*x;


%--------------------------------------
%Iteratively center the x,s pair
centered = false;
iter     = 0;

fprintf('Will use newton to find an approximately centered point \n');
while ~centered
    %Evaluate the hessian
    H = eval_socp_hessian(problem,x);
    %XXX: we use this for the step 
    %Eval the hessian at s 
    Hs = eval_socp_hessian(problem,s);

    %Evaluate the gradient
    g = eval_socp_gradient(problem,x);
    %Build the matrix
    K5= [[sparse(m,m) , A            , -b                 ,sparse(m,n)  ,sparse(m,1)];...
         [-A'         , sparse(n,n)  , c                  ,-speye(n)    ,sparse(n,1)];...
         [b'          , -c'          , sparse(1,1)        ,sparse(1,n)  ,-1         ];...
         [sparse(n,m) , sparse(mu*H) , sparse(n,1)        ,speye(n,n)   ,sparse(n,1)];...
         [sparse(1,m) , sparse(1,n)  , mu*tau^-2          ,sparse(1,n)  ,1         ]];

    %build the rhs
    rhs = [rd-(-tau*c+A'*y+s);...
           rp-(tau*b-A*x);... 
           rg-(kappa-b'*y+c'*x);...
           -s-mu*g;...
           -kappa+mu*1/tau];

    %Solve the system
    d = K5\rhs;
    dy     = d(1:m);
    dx     = d(m+1:n+m);
    dtau   = d(m+n+1);
    ds     = d(m+n+2:m+2*n+1);
    dkappa = d(m+2*n+2);

    %Barrier val
    f = eval_socp_barrier(problem,x) -log(tau)-log(kappa);

    %Calculate the H norm of dx
    x_H_norm = sqrt(dx'*H*dx);
     
    %Calculate the step size 
    lambda = sqrt(dx'*H*dx+ds'*Hs*ds+(dkappa./kappa)^2+(dtau./tau)^2);
    a      = 1/(1+lambda);
    
    if(mod(iter,10)==0||x_H_norm<1)
        fprintf('Iter %2i, lambda %3.3e, ||dx||_H^-1 %3.3e, barrier %3.3e \n',iter,lambda,x_H_norm,f); 
    end
    

    %XXXXXXXXXXXXXXXXXXXXXXXXX CONTINUE FROM HERE
    %Calculate the maximum step size 
    ratios = [-x./dx;-tau/dtau;-s./ds;-kappa/dkappa];
    max_a  = min(ratios(find(ratios>0)));
    if(max_a<a)
        fprintf('Error calculating centering step size\n');
    end

   if(x_H_norm < 1)
        centered = true;
        break;
    end
 
    x = x+a*dx;
    y = y+a*dy;
    s = s+a*ds;
    tau = tau+a*dtau;
    kappa = kappa+ a*dkappa;
    iter = iter+1;
end %End centering

%---------------------------------
%Generate the lifting directions
%Build the matrix
K5_lifting= [[sparse(m,m) , A            , -b             ,sparse(m,n)  ,sparse(m,1)];...
            [-A'         , sparse(n,n)  , c               ,-speye(n)    ,sparse(n,1)];...
            [b'          , -c'          , sparse(1,1)     ,sparse(1,n)  ,-1         ];...
            [sparse(n,m) , sparse(mu*H) , sparse(n,1)     ,-speye(n,n)  ,sparse(n,1)];...
            [sparse(1,m) , sparse(1,n)  , mu*tau^-2       ,sparse(1,n)  ,-1         ]];

%build the rhs
rhs_lifting = [rd-(-tau*c+A'*y+s);...
              rp-(tau*b-A*x);... 
              rg-(kappa-b'*y+c'*x);...
              +s+mu*g;...
              +kappa-mu*1/tau];
       
%Solve the system
d = K5_lifting\rhs_lifting;
dy_l     = d(1:m);
dx_l     = d(m+1:n+m);
dtau_l   = d(m+n+1);
ds_l     = d(m+n+2:m+2*n+1);
dkappa_l = d(m+2*n+2);

%Take the lifting step
x_lifted     = x+dx_l;
tau_lifted   = tau+dtau_l;
s_lifted     = s+ds_l;
kappa_lifted = kappa+dkappa_l;
y_lifted     = y+dy_l;

fprintf('Calculated the lifted points ||muH(x)x_l-s_l||_inf %3.3e ||mutau^-2tau_l-kappa_l||_inf %3.3e \n',...
        norm(mu*H*x_lifted-s_lifted,inf),...
        norm(mu*tau^-2*tau_lifted-kappa_lifted,inf));

%---------------------------------------------------
%Calculate the lifted step direction

%build the rhs
rhs = [rd;...
       rp;... 
       rg;...
       -s_lifted;...
       -kappa_lifted];

%Solve the system
d = K5\rhs;
dy     = d(1:m);
dx     = d(m+1:n+m);
dtau   = d(m+n+1);
ds     = d(m+n+2:m+2*n+1);
dkappa = d(m+2*n+2);

%For the LP case we can calculate the largest step to the boundary
steps      = -[x_lifted;tau_lifted;s_lifted;kappa_lifted]./[dx;dtau;ds;dkappa];
a_max_lifted  = min(steps(find(steps>0)));


%Calculate the centrality of the atd direction
lifted_gap     = zeros(sample_points,1);
lifted_points = zeros(sample_points,1);
a      = a_max_lifted;

first_found = false;
for(k=1:sample_points)
   xt  = x_lifted+a*dx;
   tt  = tau_lifted+a*dtau;
   st  = s_lifted+a*ds;
   kt  = kappa_lifted+a*dkappa;
   lifted_gap(k) = (xt'*st+kt*tt)/(nu+1);
   lifted_points(k) = a;
   a   = a*sample_bck;
   if(min(xt)>0&&min(st)>0&&min(kt,tt)>0&&~first_found)
    first_lifted = k;
    first_found = true;
   end
end




%----------------------------------------------------------
% Calculate the ATD direction

%First take an extra centering step because dx,,, etc have been calculated 
x = x+a*dx;
y = y+a*dy;
s = s+a*ds;
tau = tau+a*dtau;
kappa = kappa+ a*dkappa;

%Evaluate the hessian
H = diag(sparse((1./x).^2));
%Evaluate the gradient
g = -1./x;


fprintf('Calculated the atd points ||muH(x)x_a-s_a||_inf %3.3e ||mutau^-2tau_a-kappa_a||_inf %3.3e \n',...
        norm(mu*H*x-s,inf),...
        norm(mu*tau^-2*tau-kappa,inf));


%Build the matrix
K5= [[sparse(m,m) , A            , -b                 ,sparse(m,n)  ,sparse(m,1)];...
     [-A'         , sparse(n,n)  , c                  ,-speye(n)    ,sparse(n,1)];...
     [b'          , -c'          , sparse(1,1)        ,sparse(1,n)  ,-1         ];...
     [sparse(n,m) , sparse(mu*H) , sparse(n,1)        ,speye(n,n)   ,sparse(n,1)];...
     [sparse(1,m) , sparse(1,n)  , mu*tau^-2          ,sparse(1,n)  ,1         ]];

  
%build the rhs
rhs = [rd;...
       rp;... 
       rg;...
       -s;...
       -kappa];

%Solve the atd direction
d = K5\rhs;
dy     = d(1:m);
dx     = d(m+1:n+m);
dtau   = d(m+n+1);
ds     = d(m+n+2:m+2*n+1);
dkappa = d(m+2*n+2);

%For the LP case we can calculate the largest step to the boundary
steps      = -[x;tau;s;kappa]./[dx;dtau;ds;dkappa];
a_max_atd  = min(steps(find(steps>0)));

%Calculate the centrality of the atd direction
atd_gap    = zeros(sample_points,1);
atd_points = zeros(sample_points,1);
a      = a_max_atd;
first_found = false;
for(k=1:sample_points)
   xt  = x+a*dx;
   tt  = tau+a*dtau;
   st  = s+a*ds;
   kt  = kappa+a*dkappa;
   atd_gap(k) = (xt'*st+kt*tt)/(nu+1);
   atd_points(k) = a;
   a   = a*sample_bck;
   if(min(xt)>0&&min(st)>0&&min(kt,tt)>0&&~first_found)
    first_atd = k;
    first_found = true;
   end
end


%------------------------------------------------
% Plot the figures

figure
hold on
%Plot the gap 
plot(atd_points,atd_gap,'r');
%Mark the first found feasible point
plot(atd_points(first_atd),atd_gap(first_atd),'rd');
%Same for the lifted direction
plot(lifted_points,lifted_gap);
plot(lifted_points(first_lifted),lifted_gap(first_lifted),'bd');
hold off
title('Dual Gap vs Step Size ATD and Lifted direction');
legend('ATD','last Atd feasible','Lifted', 'last lifted feasible');
