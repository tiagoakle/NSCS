%CALL NSCS minentropy script
addpath ../matlab
% ENTROPY PROBLEM

clear; clc;
 
M = 20;
N = 50;
 
A    = randn(M,N);
xx   = 3*ones(N,1);
b    = A*xx;
d    = ones(N,1);

pars = set_default_pars_nscs_long_step();
pars.echo   = 4;
pars.beta   = 0.99;
pars.trace  = 3;
pars.secord = 1;
pars.print  = 3; 
 
[x] = minentropy_nscs(A,b,d,pars);
 
