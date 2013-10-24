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

Hx  = 2/(x'*J*x)^2*(J*x)*(J*x)' - 1/(x'*J*x)*J;
Hs  = 2/(s'*J*s)^2*(J*s)*(J*s)' - 1/(s'*J*s)*J;

H   = blkdiag(Hx,Hs);

%Test the hessian
H_eval = eval_socp_hessian(problem,[x;s]);
fprintf('||H-H_eval||: %g\n',norm(H-H_eval));

%Calculate the centering points for x,s and s,x
Jx = J*x;
xJx = x'*Jx;
Js = J*s;
sJs = s'*Js;
sx  = s'*x;
w_gamma = sqrt((1+(sx/sqrt(sJs*xJx)))/2);
w      = sqrt(sqrt((xJx/sJs)))*(1/(2*w_gamma))*(x/sqrt(xJx)+1/sqrt(sJs)*J*s);
w2     = sqrt(sqrt((sJs/xJx)))*(1/(2*w_gamma))*(J*x/sqrt(xJx)+1/sqrt(sJs)*s);

w_func = calculate_nt_scaling(problem,[x;s],[s;x]);
H_func = eval_socp_hessian(problem,w_func);

Jw    = J*w;
wJw   = w'*Jw;
Jw2   = J*w2;
w2Jw2 = w2'*Jw2; 
Hw    = 2*Jw*Jw'/(wJw)^2-1/(wJw)*J;
Hw2  = 2*Jw2*Jw2'/(w2Jw2)^2-1/(w2Jw2)*J;

fprintf('||H[x;s]-H_eval[x;s]||: %g\n',norm(H_func*[x;s]-[Hw*x;Hw2*s]));
fprintf('||H_eval[x;s]-[s;x]||: %g\n',norm(H_func*[x;s]-[s;x]));
