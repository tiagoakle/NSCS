
n = 100;
%-------------------
%Choose an initial x and s
x = randn(n(1),1);
x(1) = norm(x(2:end))+rand(1);

%Generate a feasible dual point
s    = randn(n(1),1);
s(1) = norm(s(2:end))+rand(1);

tau = rand(1);
kappa = rand(1);

%----------------------
%Choose random directions
dx     = randn(n(1),1);
ds     = randn(n(1),1);
dkappa = randn(1);
dtau   = randn(1);

[H,W,Wi,lambda] = eval_scaled_variables(n,x,s);

a_max = max_step_size(tau,kappa,dx,ds,dtau,dkappa,W,Wi,lambda)

x_h   = x+a_max*dx;
s_h   = s+a_max*ds;
tau_h   = tau+a_max*dtau;
kappa_h   = kappa+a_max*dkappa;

fprintf('tau %g, kappa %g , det x %g, det s %g\n',tau_h,kappa_h,x_h(1)-norm(x_h(2:end)),s_h(1)-norm(s_h(2:end)));
