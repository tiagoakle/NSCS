clear all
%Generate two feasible points x,s in the SOCP cone and 
%Calculate its scaling point w, evaluate the hesssian at 
%H(w) and check that H(w)x=s, Calculate W the scaling matrix
%And lambda the scaled variable


%Choose a problem size 
n = 100;
x = zeros(n,1);
s = zeros(n,1);
max_norm = 1000; %Select a bound for x0

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

slack_x = x(1)^2-norm(x(2:end))^2;
slack_s = s(1)^2-norm(s(2:end))^2;
%X and Z are not in L_n
%Sanity check
fprintf('x(1): %g, z(1): %g, x(1)^2-||x(2:)||: %g, s(1)^2-||s(2:)||: %g\n',x(1),s(1),slack_x,slack_s);

%Define the hyperbolic Householder J
J = diag(sparse([1;-ones(n-1,1)]));
%Calculate the scaled points x_hat, s_hat
scaling_x = x'*J*x;
scaling_s = s'*J*s;
s_hat     = 1/sqrt(scaling_s)*s;
x_hat     = 1/sqrt(scaling_x)*x;

%Sanity Check
fprintf('s_hat*J*s_hat %g, x_hatJx_hat %g \n',s_hat'*J*s_hat, x_hat'*J*x_hat);

%Calculate the scaling point
gamma = 0.5*(1+x_hat'*s_hat);
gamma = sqrt(gamma);

w_hat     = (x_hat+J*s_hat);
w_hat     = 0.5/gamma*w_hat;

%Evaluate the hessian at the scaled scaling point H(w)x_hat
v     = J*w_hat;
Hx    = 2*v*(v'*x_hat)-J*x_hat; %The Hessian is given by (2*vv'-J);
fprintf('||H(w)x_hat-s_hat||: %g\n',norm(Hx-s_hat));

%Calculate the unscaled scaling point
w_scal_1 = sqrt((scaling_x/scaling_s));
w_scal = sqrt(w_scal_1);
w      = w_scal*w_hat;

%Evaluate the hessian at the scaling point H(w)x
v     = J*w;
Hx    = 1/(w'*J*w)^2*2*v*(v'*x)-1/(w'*J*w)*J*x; %The Hessian is given by 1/(wJw)^2 (2*vv'-(wJw)J);

fprintf('||H(w)x-s||: %g\n',norm(Hx-s));

%Scaling matrices 



