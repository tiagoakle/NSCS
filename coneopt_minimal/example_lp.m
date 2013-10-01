% Linear problem

clear; clc;
 
M = 20;
N = 50;
 
A    = sprandn(M,N,0.8);
xx   = rand(N,1);  %Positive vector
b    = A*xx;
c    = A'*randn(M,1)+rand(N,1); %Feasible gradient

 
pars.echo   = 4;
pars.beta   = 0.99;
pars.trace  = 3;
pars.secord = 0;

pars.n = N;
pars.m = M;
%build the cone 
K.npos = N;
K.npow = 0;
K.nexp = 0;
K.nlog = 0;
K      = getbarrpar(K);

%Initial point
v0.x = ones(N,1);

R = coneopt(A,b,c,v0,K,pars);
 
