%In this test we solve the atd direction and take a step that 
% satisfies the bound that will keep mu positive
% and x feasible.
% the problem will be constructed to be linearly infeasible 

%When lambda is small enough we change all variables and achieve linear feasibility.
format compact
 close all
%Define a feasible lp
m = 10;
n = 20;

A = randn(m-1,n);
x_0 = rand(n,1);
s_0 = rand(n,1);
%Now extend A
A = [A; zeros(1,n)]
b = A*x_0; b(m) = 1
c = A'*randn(m,1)+s_0;

%Define some initial points 
x_0 = rand(n,1);
s_0 = rand(n,1);
y_0 = randn(m,1);
t_0 = 1;
k_0 = 1;
o_0 = 1;

%move the starting point far off center
x_0(1) = 1.e-10
s_0(1) = 1.e-10

%Form the initial residuals
pr = -A*x_0+t_0*b;
dr = A'*y_0+s_0-t_0*c;
gr = c'*x_0-b'*y_0+k_0;
mu_0 = x_0'*s_0+t_0*k_0;
mu_0 = mu_0/(n+1);

G =[ [zeros(m)   A    -b  pr];...
     [-A'     zeros(n) c  dr];...
     [b'        -c'    0  gr];
     [-pr'      -dr'  -gr 0 ]];

rho = (n+sqrt(n));
sigma = (n+1)/rho;

%------------------------------------------
max_iter = 1360;
phi = @(x,s,t,k) -sum(log(x))+1/(2*(x'*s+t*k)^2)*(s.*x)'*(s.*x)-2*s'*x/(x'*s+t*k)-log(t)-log(k)+rho*log(x'*s+t*k);
hnorm = @(x,s,t,k) 1/(s'*x+t*k)^2*((s.*x)'*(s.*x)+t^2*k^2) - n
cent = @(x,s,t,k) norm([x;t].*[s;k]-(x'*s+t*k)/(n+1));
logt = @(x,s,t,k) rho*log(x'*s+t*k);
fter = @(x,s,t,k) -sum(log(x))-log(t)-log(k)
logph = @(x,s,t,k)  hnorm(x,s,t,k)+logt(x,s,t,k)

cent_hist = zeros(max_iter,1);
phi_hist = zeros(max_iter,1);
hnorm_hist = zeros(max_iter,1);
logt_hist = zeros(max_iter,1);
fter_hist = zeros(max_iter,1);
logph_hist = zeros(max_iter,1);
o_hist     = zeros(max_iter,1);
linf_hist  = zeros(max_iter,1);
tomu_hist  = zeros(max_iter,1);
xds_dsx    = zeros(max_iter,1);
%-----------------------------------------


i = 0;
k = k_0;
t = t_0;
x = x_0;
y = y_0;
s = s_0;
o = o_0;
mu = mu_0;
for iter = 1:max_iter
    %Evaluate the potential
    phi_hist(iter) = phi(x,s,t,k);
    cent_hist(iter) = cent(x,s,t,k);
    hnorm_hist(iter) = hnorm(x,s,t,k);
    logt_hist(iter)  = logt(x,s,t,k);
    fter_hist(iter) = fter(x,s,t,k);
    logph_hist(iter) = logph(x,s,t,k);
    o_hist(iter)     = o;
    tomu_hist(iter)  = o*mu_0/mu;
  
    %Check that the linear constraints are satisfied
    linf_hist(iter) = norm(G*[y;x;t;o]-[zeros(m,1);s;k;0]+[zeros(m+n+1,1);mu_0*(n+1)]);
%Solve for the direction
    %Evaluate the hessian
    muH = sigma*mu*diag(x.^(-2));
    %mkot = mu*k/t; %For the symmetric form 
    mkot = sigma*mu/t^2;
    
    %Evaluate the residuals
    pra = A*x-t*b+o*pr;
    dra = -A'*y-s+t*c+o*dr;
    gra = b'*y-c'*x + o*gr -k;
    thra= -pr'*y-dr'*x-gr*t+mu_0*(n+1);
    npra = norm(pra);
    ndra = norm(dra);
    ngra = norm(gra);
    nthra= norm(thra);
    %form the rhs
    rhs = ...
    [zeros(m,1);...
    -s+mu*sigma*x.^(-1);...
    -k+mu*sigma*1/t;...
    0];

    G(m+1:m+n,m+1:m+n) = muH;
    G(m+n+1,m+n+1)     = mkot;

    %Solve for the direction
    d   = G\rhs;
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    d   = G\(-G*d+rhs);
    dy  = d(1:m);
    dx  = d(m+1:m+n);
    dt  = d(m+n+1);
    dth = d(m+n+2);
    ds  = rhs(m+1:n+m)-muH*dx;
    dk  = rhs(m+n+1)-mkot*dt;

    %Numerical erro
    n_err = G*d-rhs;
    nn_err=norm(n_err);
    neq_lin = norm([A*dx-dt*b+pr*dth;-A'*dy-ds+dt*c+dth*dr;b'*dy-c'*dx+gr*dth-k;-pr'*dy-dr'*dx-gr*dt])

    %History for the terms with dz
    xds_dsx(iter) = rho/mu*(dx'*s + dt*k + ds'*x+dk*t);
    
    %Calculate lambda
    lambdasq = (dx./x)'*(dx./x)+dt^2/t^2;
    lambda = sqrt(lambdasq);
    lambdasq_nu = lambdasq/(n+1);
    lambda_nu   = lambda/sqrt(n+1);
    alpha = 1/(1+(1-sigma)*(1+lambda_nu)+lambdasq_nu);
    alpha = min(alpha,1/(1+lambda));

    %Infeasibility of the full dual step
    dual_inf = min(s+ds);
    fprintf('Iter %i: lambda %3.3d, mu %3.3d, ||dr||: %3.3d, ||pr||: %3.3d, ||gr||: %3.3d, ||thra|| %3.3d, s_: %3.3d, norm err: %3.3d, inf s%3.3d dsTdx %3.3g',iter,lambda,mu,ndra,npra,ngra,nthra,dual_inf,nn_err,min(s),ds'*dx+dt*dk);
    
   if(min(s)<0)
    fprintf('**')
   end
%Take the step on all 
    y = y+alpha*dy;
    x = x+alpha*dx;
    s = s+alpha*ds;
    t = t+alpha*dt;
    k = k+alpha*dk;
    o = o+alpha*dth;
%Update mu 
    mu = x'*s+t*k;
    mu = mu/(n+1);
    fprintf('\n')
  
    if(mu < 1.e-16 || min(x) <0 || t<0 ||k <0 )
        if(min(x)<0)
            fprintf('Primal is infeasible!');
        end
        if(t<0)
            fprintf('Tau infeasible');
        end
        if(k<0)
            fprintf('Kappa infeasible');
        end
        if(mu<1.e-16)
            fprintf('Mu < 1.e-16')
            fprintf('Mu %f\n',mu)
        end
        phi_hist = phi_hist(1:iter);
        cent_hist = cent_hist(1:iter);
        hnorm_hist = hnorm_hist(1:iter);
        fter_hist = fter_hist(1:iter);
        logt_hist = logt_hist(1:iter);
        logph_hist= logph_hist(1:iter);
        o_hist = o_hist(1:iter);
        linf_hist = linf_hist(1:iter);
        tomu_hist = tomu_hist(1:iter);
        xds_dsx   = xds_dsx(1:iter);
        break

    end

end
plot(phi_hist)
figure 
plot(cent_hist)
title('Centrality')
figure
plot(hnorm_hist)
title('H norm')
figure
plot(logt_hist)
title('Log term')
figure
plot(fter_hist)
title('fter')
figure
plot(logph_hist)
title('\rho \log x^Ts + 1/\mu||s+\mu g(x)||')
figure
plot(o_hist)
title('theta')
figure
plot(linf_hist)
title('linear inf')
figure
plot(tomu_hist)
title('\frac{\theta\mu_0}{\mu}')
figure
plot(xds_dsx)
title('\delta x^Ts + \delta s^T x')
fprintf('\n')

