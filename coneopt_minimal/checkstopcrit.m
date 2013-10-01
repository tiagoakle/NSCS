function R = checkstopcrit(v,pars,R)

if R.stop
    % in this case, something outside says
    % that we should stop, so assign output
    % and exit with best found point so far:
    R.sol.v    = v;
    R.sol.xopt = v.x / v.tau;
    R.sol.sopt = v.s / v.tau;
    R.sol.yopt = v.y / v.tau;
    v.k = v.k + 1;
    R   = gettrace(v,pars,R);
    return;
end


if v.rPrel < pars.rhoP
    if v.rDrel < pars.rhoD
        if v.rArel < pars.rhoA

            R.stop     = true;
            R.status   = 'optimal';
            R.sol.v    = v;
            R.sol.xopt = v.x / v.tau;
            R.sol.sopt = v.s / v.tau;
            R.sol.yopt = v.y / v.tau;

            v.k = v.k + 1;
            R   = gettrace(v,pars,R);
        else
            if v.rGrel < pars.rhoG
                if v.tau < pars.rhoI*max(1,v.kappa)
                    % in this case, problem is either
                    % primal or dual infeasible
                    R.stop = true;
                    
                    %init these variables:
                    pinfeas = false; dinfeas = false;
                    
                    bty = pars.b'*v.y;
                    if  bty > 0
                        % primal infeasible
                        R.status = 'P-infeasible';
                        pinfeas  = true;
                    end
                    ctx = pars.c'*v.x;
                    if ctx < 0
                        % dual infeasible
                        R.status = 'D-infeasible';
                        dinfeas  = true;
                    end
                    R.sol.v    = v;
                    R.sol.xopt = [];
                    R.sol.sopt = [];
                    R.sol.yopt = [];
                    
                    
                    if pinfeas && dinfeas
                        tmp = abs(bty) - abs(ctx);
                        if tmp >= 0
                            R.status = 'P-infeasible';
                        else
                            R.status = 'D-infeasible';
                        end
                    end
                    if pinfeas || dinfeas
                        return;
                    end
                end
            end
        end
    end
end

% check for ill-posedness:
if v.mu < (pars.rhoM*v.mu0)
    if v.tau < pars.rhoI*min(1,v.kappa)
        R.stop     = true;
        R.status   = 'ill-posed';
        R.sol.v    = v;
        R.sol.xopt = [];
        R.sol.sopt = [];
        R.sol.yopt = [];
    end
end

% if reached maxit:
if v.k == pars.outermaxit

    R.stop   = true;
    R.status = 'maxit';

    % probably should not output anything here since
    % we do know if found sol is optimal
    R.sol.v    = v;
    R.sol.xopt = v.x / v.tau;
    R.sol.sopt = v.s / v.tau;
    R.sol.yopt = v.y / v.tau;

end


