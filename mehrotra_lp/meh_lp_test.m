clear all
close all

%This is a Mehrotra predictor-corrector implementation for LP

%SOCP solver for one cone
n = 100; %Size of the cone
m = 50; %Number of costraints
nu = n; %Complexity of the barrier

%-------
%Generate a feasible primal point
x = rand(n,1);

%Generate a feasible dual point
s = rand(n,1);

%Generate a feasible and bounded primal dual problem
A = sprandn(m,n,0.5);
b = A*x;
s = rand(n,1);
c = A'*ones(m,1) + s;
fprintf('Generated a random feasible LP with %i constraints and %i variables\n',m,n);
clear x s

%-----------------------------------
%Initialization strategy from CVXOPT
K2  = [[speye(n), A'];[A,sparse(m,m)]];
sol = K2\[zeros(n,1);b];
x   = sol(1:n);

K2  = [[speye(n), A'];[A,sparse(m,m)]];
sol = K2\[c;zeros(m,1)];
y   = sol(n+1:n+m);
s   = sol(1:n);

%%XXX:RANDOM INIT
%x   = ones(n,1);
%s   = ones(n,1);
%y   = zeros(m,1);

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
max_iter = 40;

%--------------------------------------
%Main SOCP solver iteration
for iter=1:max_iter
   
    %Evaluate the hessian and scaled variables
    %Evaluate the scaling point 
    H = diag(sparse(s./x));

    %Build the rhs term 
    %Solve the affine scalin g direction
    [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-s,-kappa,[]);
    
    %fprintf('Nres %g %g %g %g %g \n',...
    %norm(A'*dy+ds-dtau*c-rd),...
    %norm(A*dx-dtau*b),...
    %norm(b'*dy-c'*dx-dkappa-rg),...
    %norm(H*dx+ds+s),...
    %norm(kappa/tau*dtau+dkappa+kappa));
    
   % %XXX:Long step primal form of the equations?
   % H = diag(sparse(mu*(1./x).^2));
   % [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-mu*1./x,-kappa,[]);
    

    ratios = [-x./dx;-s./ds;-tau/dtau;-kappa/dkappa;1];
    a_max  = min(ratios(find(ratios>0)));

    %Now calculate the centering direction
    %-------------------------------------
    sigma = (1-a_max)^3;

%    fprintf('Sigmu %g sigma %g\n',sigmu,sigma);
%    sigma = sigmu;

    %Build the second order correction term, this should have a better form
    so       = dx.*ds./x;
    kt_so    = dtau*dkappa/tau;
      
    %so    = 0;
    %kt_so = 0;
    
    %Solve the affine scaling direction
    [dy_c,dx_c,dtau_c,ds_c,dkappa_c,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,...
                                                                                            (1-sigma)*rp,...
                                                                                            (1-sigma)*rd,...
                                                                                            (1-sigma)*rg,...
                                                                                            -s    +sigma*mu*1./x  -so,...
                                                                                            -kappa+sigma*mu*1/tau -kt_so,...
                                                                                            slv_aug);


    %Decrement 
    ratios = [-x./dx_c;-s./ds_c;-tau/dtau_c;-kappa/dkappa_c;1];
    a_max  = min(ratios(find(ratios>0)));
    a       = a_max*0.98; 
    % a       = min(a,1); 
    
    x_old = x;
    y_old = y;
    s_old = s;
    t_old = tau;
    k_old = kappa;

    y         = y+a*dy_c;
    x         = x+a*dx_c;
    s         = s+a*ds_c;
    tau       = tau+a*dtau_c;
    kappa     = kappa+a*dkappa_c;

    mu        = x'*s+tau*kappa;
    mu        = mu/(nu+1);
    
    %Feasiblity check
    x_slack = min(x);
    s_slack = min(s);
    if(x_slack<0||s_slack<0||tau<0||kappa<0)
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



%now call ecos 
G = -speye(n);
h = zeros(n,1);
dims.l = n;
dims.q = [];
[x_ecos,y_ecos,info_ecos,s_ecos,z_ecos] = ecos(c,G,h,dims,A,b);
