
%Define the problem size
m = 5;
n_1 = 10;
n_2 = 20;

%-------
%Generate a feasible primal point
x = randn(n_1,1);
x(1) = norm(x(2:end))+rand(1);
x_feas_1 = x;

x = randn(n_2,1);
x(1) = norm(x(2:end))+rand(1);
x_feas_2 = x;


%Generate a feasible dual point
s = randn(n_1,1);
s(1) = norm(s(2:end))+rand(1);
s_feas_1 = s;

s = randn(n_2,1);
s(1) = norm(s(2:end))+rand(1);
s_feas_2 = s;


%Generate a feasible and bounded primal dual problem
A = sprandn(m,n_1+n_2,0.5);
b = A*[x_feas_1;x_feas_2];
s = rand(n,1);
c = A'*ones(m,1) + [s_feas_1;s_feas_2];
fprintf('Generated a random feasible SOCP with %i constraints and %i variables\n',m,n);

socp_cones = [n_1,n_2];
[x,y,s,info] = socp_solver(A,b,c,socp_cones);
