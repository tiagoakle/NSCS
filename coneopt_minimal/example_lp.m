clear all
clc

%Build an lp problem
m = 50;
n = 100;
nf = 0;
nc = n-nf;

A  = randn(m,n);
x  = [randn(nf,1);rand(nc,1)];
b  = A*x;

c          = A'*randn(m,1);
c(nf+1:n)  = c(nf+1:n) + rand(n-nf,1);

%Set parameters
pars.n = n;
pars.m = m;
pars.cnbfgsstps = 0;
pars.echo = 4;
v0.x     = rand(nc,1);

% build cone:
K.npos = nc;
K.npow = 0;
K.nexp = 0;
K.nlog = 0;
K      = getbarrpar(K);

% call to coneopt:
R = coneopt(A,b,c,v0,K,pars);


