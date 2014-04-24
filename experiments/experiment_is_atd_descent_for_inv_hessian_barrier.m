 %Experiment is ATD a potential reduction direction for the 
 %function defined with no dual barrier but the term 1/2mu * s'H(x)^{-1}s - 2s'x instead
% Generate a random lp and calculate the ATD direction and the symmetric direction
% plot the value of Psi(x,s) = \rho \log(x's) + f(x) + 1/2mu * s'H(x)^{-1}s - 2s'x. Allong the atd and the symmetric direction
% Rho should be chosen such that sigma = nu/rho
close all;
clear all;

%Problem size
n = 100;
m = 10;
nu = n;
sigma = 0.1;
rho = nu/sigma;

%Graph parameters
samples = 1000;

%Generate a random x,s 
x = rand(n,1);
s = rand(n,1);

%Make an entry particulary off center
s(1) = 1e-10;
x(1) = 1e-10;

mu = x'*s/nu;

%Generate an LP for which these are feasible
A = randn(m,n);
y = randn(m,1);
b = A*x;
c = A'*y+s;

%------------------------------------------
%Calculate the ATD direction

%Evaluate the hessiand and gradient
H = diag((1./x).^2);
g = -(1./x);

%Calculate the ATD direction 
K = [[-mu*H A'];[A zeros(m,m)]];
rhs = [s+sigma*mu*g;zeros(m,1)];
d = K\rhs;
dx_atd = d(1:n);
dy_atd = d(n+1:n+m);
ds_atd = -A'*dy_atd;

%-----------------------------------------
%Calculate the symmetric direction 

%Evaluate the hessian and gradient
H = diag(s./x);
g = -(1./x);

%Calculate the ATD direction 
K = [[-mu*H A'];[A zeros(m,m)]];
rhs = [s+sigma*mu*g;zeros(m,1)];
d = K\rhs;
dx_sym = d(1:n);
dy_sym = d(n+1:n+m);
ds_sym = -A'*dy_sym;


%Find the largest step admissible for both directions
ratios = [-x./dx_sym;-s./ds_sym;-x./dx_atd;-s./ds_atd];
%ratios = [-x./dx_sym;-x./dx_atd];
max_step = min(ratios(find(ratios>0)));
if(isempty(max_step))
    max_step_atd = 1;
end


%---------------------------------------
%Generate the plots

%Store the results here
alphas   = zeros(samples,1);
vals_atd = zeros(samples,1);
vals_sym = zeros(samples,1);

%Define the potential function
psi = @(u,v)(rho*log(u'*v)-sum(log(u)) + 0.5/(u'*v)*(v'*diag(u.^2)*v)-2*v'*u );

alphas = max_step*[0:samples-1]./(samples-1)*0.99995;
for j = 1:samples
    vals_atd(j) = psi(x+alphas(j)*dx_atd,s+alphas(j)*ds_atd);
    vals_sym(j) = psi(x+alphas(j)*dx_sym,s+alphas(j)*ds_sym);
end

figure
hold on
plot(alphas,vals_atd,'r');
%plot(alphas,vals_sym,'g');
hold off;


%Calculate the inner product with the gradient 
grad_psi = [rho/mu*s-1./x;rho/mu*x-1./s];
%gpsiTatd = grad_psi'*[dx_atd;ds_atd]
%gpsiTsym = grad_psi'*[dx_sym;ds_sym]
