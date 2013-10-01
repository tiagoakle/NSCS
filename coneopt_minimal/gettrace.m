function R = gettrace(v,pars,R)

if pars.trace > 0
    R.trace.tau(v.k)         = v.tau;
    R.trace.kappa(v.k)       = v.kappa;
    R.trace.mu(v.k)          = v.mu;
    R.trace.a(v.k)           = v.a;
    R.trace.rPrel(v.k)       = v.rPrel;
    R.trace.rDrel(v.k)       = v.rDrel;
    R.trace.rGrel(v.k)       = v.rGrel;
    R.trace.rArel(v.k)       = v.rArel;
    R.trace.ncentstepse(v.k) = v.i;
    R.trace.lam(v.k)         = v.lam;
    R.trace.feas(v.k)        = v.feas;
end
if pars.trace > 1
    R.trace.x(:,v.k) = v.x;
    R.trace.s(:,v.k) = v.s;
    R.trace.y(:,v.k) = v.y;
end