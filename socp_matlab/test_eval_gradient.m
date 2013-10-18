clear all
%Test the eval gradient method

%Choose a problem size 
n = 10;
x = zeros(n,1);
s = zeros(n,1);
max_norm = 10; %Select a bound for x0

%Choose a random x,s
%Pick a pair of random points inside the unit sphere of 
%size one

d = randn(n-1,1); %Pick a direction
d = d/norm(d);
x(1)   = max_norm*rand(1);
x(2:n) = x(1)*d*sqrt(rand(1)); % Pick a distance from the origin

d = randn(n-1,1); %Pick a direction
d = d/norm(d);
s(1)   = max_norm*rand(1);
s(2:n) = s(1)*d*sqrt(rand(1)); % Pick a distance from the origin

J = diag(sparse([1;-ones(n-1,1)]));


%Calculate the gradients at x and s
gx = -1/(x'*J*x)*J*x;
gs = -1/(s'*J*s)*J*s;

g  = [gx;gs];

%Form the problem structure and call the method
problem = struct;
problem.n_soc_cones = 2;
problem.soc_cones   = [n,n];

g_eval = eval_socp_gradient(problem,[x;s]);
fprintf('||g-g_eval||: %g\n',norm(g-g_eval));
