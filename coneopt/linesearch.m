function [a,R] = linesearch(v,K,pars,R)

if pars.useamax
    amax = 1.0;
    a    = getamax(v.x,v.s,v.tau,v.kappa,v.dx,v.ds,...
        v.dtau,v.dkappa,amax,pars.lscaff,v,K,pars);
else
    % set maximal step length:
    a0 = 1.0;
    a  = a0;
    if v.dkappa < 0
        kapmax = -v.kappa/v.dkappa;
        a = min(a,kapmax);
    end
    if v.dtau < 0
        taumax = -v.tau/v.dtau;
        a = min(a,taumax);
    end
    % if either kap or tau is blocking, multiply by eta
    % so we do not hit boundary:
    if a < a0
        a = pars.eta*a;
    end
end


R.block = 'non';
if pars.trace > 0
    R.dat.amaxaff(v.k) = a;
end

% count number of bisections:
nsect = 0;

for j = 1:pars.lsmaxit
    
    % point to try:
    xa     = v.x     + a * v.dx;
    sa     = v.s     + a * v.ds;
    taua   = v.tau   + a * v.dtau;
    kappaa = v.kappa + a * v.dkappa;
    
    % new duality gap:
    dga    = xa'*sa + taua*kappaa;
    mua    = dga / (K.nu + 1);
    
    % evaluate barriers at new point:
    % check only feasibility, so want = [-1,-1,-1]:
    FP = BarrFuncP(xa,K,[1,1,1]);
    FD = BarrFuncD(sa,K,[1,-1,-1]);
    
    dosect = false;
    if FP{4} < 0
        dosect  = true;
        R.block = 'pf';
    elseif FD{4} < 0 && pars.checkDpt
        dosect  = true;
        R.block = 'df';
    else
        psi       = sa + mua*FP{2};
        centmeas1 = norm( [sa./FP{2};-taua*kappaa] + mua, inf);
        centmeas2 = norm(psi,inf)/norm(FP{2},inf);
        centmeas2 = max(centmeas2,abs(taua*kappaa-mua));
        centmeas3 = abs(xa'*psi);
        centmeas4 = norm(psi);
        centmeas5 = sqrt(psi'*(FP{3}\psi));
        if pars.centmeastype == 1
            centmeas = centmeas1;
        elseif pars.centmeastype == 2
            centmeas = centmeas2;
        elseif pars.centmeastype == 3
            centmeas = centmeas3;
        elseif pars.centmeastype == 4
            centmeas = centmeas4;
        elseif pars.centmeastype == 5
            centmeas = centmeas5;
        end
        if centmeas > mua*pars.theta
            dosect  = true;
            R.block = 'ce';
        end
    end
    
    if dosect
        a     = a*pars.lscaff;
        nsect = nsect + 1;
    else
        break;
    end
    
end

% safety:
if j == pars.lsmaxit
    xa = v.x + a * v.dx;
    FP = BarrFuncP(xa,K,[1,1,-1]);
    if FP{4} < 0
        error(['linesearch: failed to find feasible point.',...
            ' pars.lscaff too close to 1 ???']);
    end
end



