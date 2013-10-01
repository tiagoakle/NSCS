function [v,R] = slvcentdir_coneopt(v,rs,K,pars,R)

% get the current right hand side:
r1  = rs.r1;
r2  = rs.r2;
r3  = rs.r3;
r4  = rs.r4;
r5  = rs.r5;

%testing:
%ck = v.mu*(1/v.kappa);
ck = v.tau;


[d,v.CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
         r1,r2,r3,r4,r5,pars);

%[d,v.CF] = slvhomkkt(v.F{3},v.mu,pars.A,pars.b,pars.c,ck,v.kappa,...
%    r1,r2,r3,r4,r5,pars);

if d{6} > 0
    R.status = 'rounding (cent)';
    R.stop   = true;
    return;
end

% extract sol:
v.dxc     = d{1};
v.dtauc   = d{2};
v.dyc     = d{3};
v.dsc     = d{4};
v.dkappac = d{5};

% counting:
R.dat.nkktsolves = R.dat.nkktsolves  + 1;


