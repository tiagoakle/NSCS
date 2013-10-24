clear all
close all

%SOCP solver for one cone
n = 20; %Size of the cone
m = 10; %Number of costraints
n = sum(problem.soc_cones);
nu= 2*problem.n_soc_cones; %Parameter of the problem

%-------
%Generate a feasible primal point
x1 = randn(problem.soc_cones(1),1);
x2 = randn(problem.soc_cones(2),1);
x3 = randn(problem.soc_cones(3),1);
x1(1) = norm(x1(2:end))+rand(1);
x2(1) = norm(x2(2:end))+rand(1);
x3(1) = norm(x3(2:end))+rand(1);
x_feas = [x1;x2;x3];

%Generate a feasible dual point
s1 = randn(problem.soc_cones(1),1);
s2 = randn(problem.soc_cones(2),1);
s3 = randn(problem.soc_cones(3),1);
s1(1) = norm(s1(2:end))+rand(1);
s2(1) = norm(s2(2:end))+rand(1);
s3(1) = norm(s3(2:end))+rand(1);
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
x1(1) = norm(x1(2:end))+rand(1);
x2(1) = norm(x2(2:end))+rand(1);
x3(1) = norm(x3(2:end))+rand(1);
x  = [x1;x2;x3];

%Generate a feasible dual point
s1 = randn(problem.soc_cones(1),1);
s2 = randn(problem.soc_cones(2),1);
s3 = randn(problem.soc_cones(3),1);
s1(1) = norm(s1(2:end))+rand(1);
s2(1) = norm(s2(2:end))+rand(1);
s3(1) = norm(s3(2:end))+rand(1);
s  = [s1;s2;s3];
tau = 1;
kappa = 1;
fprintf('Selected random x,s \n');
y  = ones(m,1);

%------Calculate the initial mu
mu = (s'*x+tau*kappa)/(nu+1);

%--------------------
%Evaluate the residuals
rp = tau*b-A*x;
rd = -tau*c+A'*y+s;
rg = kappa-b'*y+c'*x;


%--------------------------------------
%Main SOCP solver iteration
solved = false;
iter     = 0;

while ~solved
   
    %Evaluate the hessian
    H = eval_socp_hessian(problem,x);
    
    %Evaluate the gradient
    g = eval_socp_gradient(problem,x);
  
    %Build the matrix
    K5= [[sparse(m,m) , A            , -b                 ,sparse(m,n)  ,sparse(m,1)];...
         [-A'         , sparse(n,n)  , c                  ,-speye(n)    ,sparse(n,1)];...
         [b'          , -c'          , sparse(1,1)        ,sparse(1,n)  ,-1         ];...
         [sparse(n,m) , sparse(mu*H) , sparse(n,1)        ,speye(n,n)   ,sparse(n,1)];...
         [sparse(1,m) , sparse(1,n)  , mu*tau^-2          ,sparse(1,n)  ,1         ]];

    %build the rhs for the affine direction
    rhs = [rd;...
           rp;... 
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

    %Calculate the maximum step size 
    a_max = 1; %Where to store the max steps
    i     = 1;
    for(k = 1:problem.n_soc_cones)

        if(-x(i)/dx(i)>0)
            a_max = min(a_max,-dx(i)/x(i));
        end

        b_p = 2*(dx(i)*x(i) - dx(i+1:i+problem.soc_cones(k)-1)'*x(i+1:i+problem.soc_cones(k)-1));
        c_p = x(i)^2-norm(x(i+1:i+problem.soc_cones(k)-1))^2;
        a_p = dx(i)^2-norm(dx(i+1:i+problem.soc_cones(k)-1))^2;
        deter = b_p^2-4*a_p*c_p;
        if(deter>0)
            sd = sqrt(deter);
            pr = (-b_p+sd)/(2*a_p);
            nr = (-b_p-sd)/(2*a_p);
            if(pr<0&&nr<0)
                a_max = min(a_max,inf);
            elseif(nr>0&&pr>0)
                a_max = min(a_max,min(nr,pr));
            else
                a_max = min(a_max,max(nr,pr));
            end
        end     

        %DO the same for s

        if(-s(i)/ds(i)>0)
            a_max = min(a_max,-ds(i)/s(i));
        end

        b_p = 2*(ds(i)*s(i) - ds(i+1:i+problem.soc_cones(k)-1)'*s(i+1:i+problem.soc_cones(k)-1));
        c_p = s(i)^2-norm(s(i+1:i+problem.soc_cones(k)-1))^2;
        a_p = ds(i)^2-norm(ds(i+1:i+problem.soc_cones(k)-1))^2;
        deter = b_p^2-4*a_p*c_p;
        
        if(deter>0)
            sd = sqrt(deter);
            pr = (-b_p+sd)/(2*a_p);
            nr = (-b_p-sd)/(2*a_p);
            if(pr<0&&nr<0)
                a_max = min(a_max,inf);
            elseif(nr>0&&pr>0)
                a_max = min(a_max,min(nr,pr));
            else
                a_max = min(a_max,max(nr,pr));
            end
        end     
      
        i = i+problem.soc_cones(k);
    end
     
    ratios = [-tau/dtau,-kappa/dkappa,a_max];
    a_max  = min(ratios(find(ratios>0)));

    %Now calculate the centering and corrector direction
    %--------------------------------------------------
    sigma = (1-a_max)^3;

     %build the rhs for the centering direction
    rhs = [rd-(-tau*c+A'*y+s);...
           rp-(tau*b-A*x);... 
           rg-(kappa-b'*y+c'*x);...
           -s-sigma*mu*g;...
           -kappa+sigma*mu/tau];
 
   
    %Solve the system
    d = K5\rhs;
    dy_c     = d(1:m);
    dx_c     = d(m+1:n+m);
    dtau_c   = d(m+n+1);
    ds_c     = d(m+n+2:m+2*n+1);
    dkappa_c = d(m+2*n+2);

    %dx       = dx+dx_c;
    %dy       = dy+dy_c;
    %ds       = ds+ds_c;
    %dtau     = dtau+dtau_c;
    %dkappa   = dkappa+dkappa_c;
    
    dx       = dx_c;
    dy       = dy_c;
    ds       = ds_c;
    dtau     = dtau_c;
    dkappa   = dkappa_c;


    %Calculate the maximum step size 
    a_max = 1; %Where to store the max steps
    i     = 1;
    for(k = 1:problem.n_soc_cones)

        if(-x(i)/dx(i)>0)
            a_max = min(a_max,-dx(i)/x(i));
        end

        b_p = 2*(dx(i)*x(i) - dx(i+1:i+problem.soc_cones(k)-1)'*x(i+1:i+problem.soc_cones(k)-1));
        c_p = x(i)^2-norm(x(i+1:i+problem.soc_cones(k)-1))^2;
        a_p = dx(i)^2-norm(dx(i+1:i+problem.soc_cones(k)-1))^2;
        deter = b_p^2-4*a_p*c_p;
        if(deter>0)
            sd = sqrt(deter);
            pr = (-b_p+sd)/(2*a_p);
            nr = (-b_p-sd)/(2*a_p);
            if(pr<0&&nr<0)
                a_max = min(a_max,inf);
            elseif(nr>0&&pr>0)
                a_max = min(a_max,min(nr,pr));
            else
                a_max = min(a_max,max(nr,pr));
            end
        end     

        %DO the same for s

        if(-s(i)/ds(i)>0)
            a_max = min(a_max,-ds(i)/s(i));
        end

        b_p = 2*(ds(i)*s(i) - ds(i+1:i+problem.soc_cones(k)-1)'*s(i+1:i+problem.soc_cones(k)-1));
        c_p = s(i)^2-norm(s(i+1:i+problem.soc_cones(k)-1))^2;
        a_p = ds(i)^2-norm(ds(i+1:i+problem.soc_cones(k)-1))^2;
        deter = b_p^2-4*a_p*c_p;
        
        if(deter>0)
            sd = sqrt(deter);
            pr = (-b_p+sd)/(2*a_p);
            nr = (-b_p-sd)/(2*a_p);
            if(pr<0&&nr<0)
                a_max = min(a_max,inf);
            elseif(nr>0&&pr>0)
                a_max = min(a_max,min(nr,pr));
            else
                a_max = min(a_max,max(nr,pr));
            end
        end     
      
        i = i+problem.soc_cones(k);
    end
     
    ratios = [-tau/dtau,-kappa/dkappa,a_max];
    a_max  = min(ratios(find(ratios>0)));
    if(isempty(a_max)) a_max = 1; end
    a         = 0.995*a_max;

    y         = y+a*dy;
    x         = x+a*dx;
    s         = s+a*ds;
    tau       = tau+a*dtau;
    kappa     = kappa+a*dkappa;
 
    mu        = x'*s+tau*kappa;
    mu        = mu/(nu+1);


    %Feasiblity check
    i     = 1;
    for(k = 1:problem.n_soc_cones)
        x_slack = x(i)^2-norm(x(i+1:i+problem.soc_cones(k)-1))^2;
        s_slack = s(i)^2-norm(s(i+1:i+problem.soc_cones(k)-1))^2;
        if(x_slack<0||s_slack<0||x(i)<0||s(i)<0||tau<0||kappa<0)
            fprintf('Infeasible centering step \n');
            return;
        end
        i = i+problem.soc_cones(k);
    end



    %Evaluate the residuals
    rp = tau*b-A*x;
    rd = -tau*c+A'*y+s;
    rg = kappa-b'*y+c'*x;

    fprintf('%2i a %3.3g pr %3.3g dr %3.3g gr %3.3g mu %3.3g \n',iter,a,norm(rp),norm(rd),norm(rg),mu);

    iter = iter+1;
end %End centering


