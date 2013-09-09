function [v,R] = slvatd_coneopt(v,K,pars,R)

% solve for
% approximate tangent direction

r4 = -v.kappa*v.tau;
r5 = -v.s;

[d,CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
    -v.rP,-v.rD,-v.rG,r4,r5,pars);

% counting:
R.dat.nkktsolves = R.dat.nkktsolves  + 1;

if d{6} > 0
    R.status = 'rounding (atd)';
    R.stop   = true;
    return;
end

v.dx     = d{1};
v.dtau   = d{2};
v.dy     = d{3};
v.ds     = d{4};
v.dkappa = d{5};

v.dtauaff   = v.dtau;
v.dkappaaff = v.dkappa;

end






