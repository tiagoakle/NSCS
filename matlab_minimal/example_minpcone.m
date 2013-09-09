% P-CONE PROBLEM
clear; clc;

M = 20;
N = 70;
A = 2*randn(M,N);
b = randn(M,1);
p = 5/3;


% set some parameters
pars.beta     = 0.99;
pars.beta2    = 0.5;
pars.scalpo   = 0;
pars.echo     = 4;
pars.permuteM = 0;
pars.theta    = 0.95;
pars.useamax  = 0;
pars.bfgsstop = 0;
pars.rhoP     = 1e-6;
pars.cnbfgsstps = 0;
pars.secord2    = 0;

[x1,f1,R1] = minpcone(A,b,p,pars);
