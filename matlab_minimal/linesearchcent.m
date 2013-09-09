function v = linesearchcent(v,K,pars)

if pars.useamax
    amax = 1.0;
    am   = getamax(v.x,v.s,v.tau,v.kappa,v.dxc,v.dsc,...
        v.dtauc,v.dkappac,amax,pars.lsccent,v,K,pars);
    a    = am;
else
    % set maximal step length:
    a0 = pars.eta*1.0;
    a  = a0;
    if v.dkappac < 0
        kapmax = -v.kappa/v.dkappac;
        a = min(a,kapmax);
    end
    if v.dtauc < 0
        taumax = -v.tau/v.dtauc;
        a = min(a,taumax);
    end
    % if either kap or tau is blocking, multiply by eta
    % so we do not hit boundary:
    if a < a0
        a = pars.eta*a;
    end

end

% init:
objvalb = v.centres;
nsect   = 0;
block   = 'non';

for j = 1:pars.lsmaxit

    % point to try:
    xa     = v.x     + a * v.dxc;
    ya     = v.y     + a * v.dyc;
    sa     = v.s     + a * v.dsc;
    taua   = v.tau   + a * v.dtauc;
    kappaa = v.kappa + a * v.dkappac;

    % objective value at xa:
    FPa = BarrFuncP(xa,K,[1,1,1]);

    dosect = false;
    if FPa{4} < 0  %Check for infeasibility
        dosect = true;
        block  = 'pfeas';
    elseif pars.dualcent && FDa{4} < 0
        dosect = true;
        block  = 'dfeas';
    else

        %r2      = v.rD0;
        r5      = -(-pars.A'*ya + v.mu*FPa{2} + taua*pars.c);
        r2cent  = -sa - v.mu*FPa{2};
        tmp     = r2cent/v.mu;
        [L,pp]  = chol(FPa{3});
        if pp > 0
            objvalt = 99;
        else
            tmp     = L'\tmp;
            objvalt = sqrt(tmp'*tmp);
        end

        if objvalt < objvalb
            objvalb = objvalt; % new best objval
        else
            dosect = true;
            block  = 'objval';
        end
    end

    if dosect
        a     = a * pars.lsccent;
        nsect = nsect + 1;
    else
        break;
    end

    % safety 1:
    if a < pars.stpszmin
        error('linesearchcent: step size smaller than minimum');
        break;
    end

    % safety 2:
    if nsect > pars.nmaxsect 
        if FPa{4} < 0
            error('linesearchcent: failed to find feasible point');
        end
        break;
    end

end

%if a < 1e-6
%    fprintf([block,'\n']);
%end
v.ac         = a;
v.centobjval = objvalb;
v.r5cent     = r5;


