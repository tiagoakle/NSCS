clear all
close all
%The objective of this experiment is to compare the rate at which 
% x'*s decays along two different directions:
% The atd direction for (x,s) calculated from muH(x)dx+dx=-s
% The centered direction for (x+dxc,s+dxc) calculated from muH(x)dx+ds=-s-dxl
% where muH(x)dxl+dsl=s+mug(x) and ||dxl||_H(x) < 1

% The experiment consists of the following: 
% 1- Generate a random LP
% 2- Generate a cone feasible x0,tau0
% 3- Generate a cone feasible s0,kappa0
% 4- Set y0 = 0
% 5-Calculate the residuals rp = taub-Ax, rd = -tauc+s+A'y, rg = kappa -b'y+c'x

% 6- Repeatedly solve the equations
% G{dyc\\dxc} - {0\\dsc} = {rd,rp,rg}
% \muH(x)dxc + dsc      = -s-mug
% and take a step x=x+adx, s=s+ads
% until ||dxc||_H < mu

% 7- Calculate dx,ds for atd as described above
% 8- Calculate the centered points x-dxc s+dsc, and dx,ds for the centered direction
% 9- Plot the results

%----------
%Experiment parameters
m = 10;
n = 20;
nu=n; %Parameter of the cone

%Plotting parameters
sample_points = 1000;
sample_bck    = 0.99; %The samples are collected at the points (sample_bck)^k

%-------
%Generate a random strictly feasible LP
A = randn(m,n);
b = A*rand(n,1);
s = rand(n,1);
c = A'*ones(m,1) + rand(n,1);
fprintf('Generated a random feasible LP with %i constraints and %i variables\n',m,n);

%-------------------
%Choose an initial x and s
s= rand(n,1);
kappa = rand(1,1);
x = rand(n,1);
tau = rand(1,1);
mu = 0.05*(s'*x+tau*kappa)/(nu+1);
y  = ones(m,1);
fprintf('Selected random x,s \n');

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
    H = diag(sparse((1./x).^2));
    %Evaluate the gradient
    g = -1./x;
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
    f = sum(-log(x)-log(s)-log(tau)-log(kappa));

    %Calculate the H norm of dx
    x_H_norm = sqrt(dx'*H*dx);

    %Calculate the step size 
    lambda = sqrt(sum((dx./x).^2+(ds./s).^2)+(dkappa./kappa)^2+(dtau./tau)^2);
    a      = 1/(1+lambda);
    
    if(mod(iter,10)==0||x_H_norm<1)
        fprintf('Iter %2i, lambda %3.3e, ||dx||_H^-1 %3.3e, barrier %3.3e \n',iter,lambda,x_H_norm,f); 
    end
    
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
