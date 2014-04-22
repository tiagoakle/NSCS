%In this test we solve the atd direction
%With a term for the dual residual and change all variables but s,
%When lambda is small enough we change all variables and achieve linear feasibility.



%Define a feasible lp
m = 20;
n = 300;

A = randn(m,n);
x_0 = rand(n,1);
s_0 = rand(n,1);
b = A*x_0;
c = A'*randn(m,1)+s_0;

%Define some initial points 
x_0 = rand(n,1);
s_0 = rand(n,1);
y_0 = randn(m,1);
t_0 = 1;
k_0 = 1;
o_0 = 1;

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
sig = n/rho;

%------------------------------------------
max_iter = 1360;
phi = @(x,s,t,k) -sum(log(x))-sum(log(s))-log(t)-log(k)+rho*log(x'*s+t*k);
phi_hist = zeros(max_iter,1);
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
%Solve for the direction
    %Evaluate the hessian
    muH = mu*diag(x.^(-2));
    mkot = mu*k/t; 
    %Evaluate the potential
    phi_hist(iter) = phi(x,s,t,k);
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
    -dra-s+mu*sigma*x.^(-1);...
    -k+mu*sigma*1/t;...
    0];

    G(m+1:m+n,m+1:m+n) = muH;
    G(m+n+1,m+n+1)     = mkot;

    %Solve for the direction
    d   = G\rhs;
    dy  = d(1:m);
    dx  = d(m+1:m+n);
    dt  = d(m+n+1);
    dth = d(m+n+2);
    ds  = rhs(m+1:n+m)+dra-muH*dx;
    dk  = rhs(m+n+1)-mkot*dt;

    %Numerical erro
    n_err = G*d-rhs;
    nn_err=norm(n_err);

%Calculate lambda
    lambda = (dx./x)'*(dx./x)+dt^2*k/t+dk^2*t/k;
    lambda = sqrt(lambda);
    alpha = 1/(1+lambda);

    %Infeasibility of the full dual step
    dual_inf = min(s+ds);
    fprintf('Iter %i: lambda %3.3d, mu %3.3d, ||dr||: %3.3d, ||pr||: %3.3d, ||gr||: %3.3d, ||thra|| %3.3d, s_: %3.3d, norm err: %3.3d',iter,lambda,mu,ndra,npra,ngra,nthra,dual_inf,nn_err);

% If lambda < sigma take a full step on all variables 
    if lambda < sigma

        alpha = 1; %overrides the previous value 
        s = s+ds;

        derr = -A'*dy-ds+dt*c+dth*dr+dra;
        nderr = norm(derr)
        fprintf(' ** ');
    end
    
%Take the step on all but s
    y = y+alpha*dy;
    x = x+alpha*dx;
    t = t+alpha*dt;
    k = k+alpha*dk;
    o = o+alpha*dth;
%Update mu 
    mu = x'*s+t*k;
    mu = mu/(n+1);
    fprintf('\n')

    if(mu < 1.e-10 || min(x) <0 || t<0 ||k <0|| min(s) <0)
        if(min(x)<0)
            fprintf('Primal is infeasible!');
        end
        if(t<0)
            fprintf('Tau infeasible');
        end
        if(k<0)
            fprintf('Kappa infeasible');
        end
        if(min(s)<0)
            fprintf('s infeasible');
        end
        phi_hist = phi_hist(1:iter);
        break

    end
end

fprintf('\n')
plot(phi_hist)
