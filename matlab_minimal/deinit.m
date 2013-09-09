function [v,K,pars,R]  = deinit(v,K,pars,R)

if isfield(R,'sol')
    R.dat.avgcentsteps = R.dat.ncentsteps / R.sol.v.k;
else
    R.dat.avgcentsteps = inf;
end
R.dat.tt           = toc-R.dat.t1;
[R.m,R.n]          = size(pars.A);