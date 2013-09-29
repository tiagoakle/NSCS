%Loads the problem in ./test_data/minentropy/matlab_copy
% and runs it in coneoptl
addpath '../interface'
addpath '.'
load('../C/test/test_data/minentropy/matlab_copy');
%Remove the permutation in A
A = sparse(m,n,nnzA);
A(:,permute) = AA;
c(permute)   = c;
% build cone:
K.npos = 0;
K.npow = 0;
K.nexp = N;
K.nlog = 0;
K      = getbarrpar(K);
 
pars.echo   = 4;
pars.beta   = 0.99;
pars.trace  = 3;
pars.secord = 0;

pars.n = n;
pars.m = m;
v0.x  = [u0;v00;x0];
% call to coneopt:
R = coneopt(A,bb,c,v0,K,pars);


