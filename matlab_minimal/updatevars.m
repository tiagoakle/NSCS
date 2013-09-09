function v = updatevars(v,K)

% store old:
v.xprev     = v.x;
v.gprev     = v.F{2};
v.tauprev   = v.tau;
v.kappaprev = v.kappa;

% take step:
v.x     = v.x     + v.a * v.dx;
v.tau   = v.tau   + v.a * v.dtau;
v.y     = v.y     + v.a * v.dy;
v.s     = v.s     + v.a * v.ds;
v.kappa = v.kappa + v.a * v.dkappa;

% update other quantities:
v.dgap  = v.x'*v.s + v.tau*v.kappa;
v.mu    = v.dgap / (K.nu + 1);

% "feas" measure, see Sturm:
v.feas  = v.dtauaff/v.tauprev - v.dkappaaff/v.kappaprev;