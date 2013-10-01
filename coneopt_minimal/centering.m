function [v,R] = centering(v,K,pars,R)

v.i     = 0;
v.lam   = 99;
v.jt    = 0;


fullcentsteps = 0;
tbfgssteps    = 0;
stopcentering = false;

if ~R.stop
    
%    v.rD0  = -pars.A'*v.y - v.s + v.tau*pars.c;
    
    for i = 1:pars.innermaxit
        v.i = i;
        
        % evalute barrier at new x:
        v.F = BarrFuncP(v.x,K,[1,1,1]);
       
        % right-hand sides:
        rs.r1 = zeros(size(pars.A,1),1);
        rs.r2 = zeros(size(pars.A,2),1);
        rs.r3 = 0;
        rs.r4 = v.mu   - (v.tau*v.kappa);
        rs.r5 = - v.s-v.mu*v.F{2}; %        - (-pars.A'*v.y + v.mu*v.F{2} + v.tau*pars.c);
 

        % compute first centering measure:
        r2cent     = rs.r5;
        tmp        = r2cent/v.mu;
        %Cholesky of the hessian
        [L,pp]     = chol(v.F{3});
        if pp > 0      % indicates v.F{3} is not pd
            v.centres  = 99;
        else
            tmp        = L'\tmp;
            v.centres  = sqrt(tmp'*tmp);
        end
         
        if stopcentering
            break;
        end
        fullcentsteps = fullcentsteps + 1;
        
        %counting full centering steps:
        R.dat.ncentsteps = R.dat.ncentsteps + 1;
        
        % solve for centering direction:
        [v,R] = slvcentdir_coneopt(v,rs,K,pars,R);
        
        % stop if rounding is blocking progress:
        if R.stop, R = checkstopcrit(v,pars,R); break; end
        
        % line search:
        v = linesearchcent(v,K,pars);
        
        % store previous step before taking new step (for bfgs)
        v.xprev = v.x;
        v.gprev = v.F{2};
        
        % take step:
        v.x     = v.x     + v.ac * v.dxc;
        v.tau   = v.tau   + v.ac * v.dtauc;
        v.y     = v.y     + v.ac * v.dyc;
        v.kappa = v.kappa + v.ac * v.dkappac;
        v.s     = v.s     + v.ac * v.dsc;
        
        % current duality gap:
        v.cdgap = v.x'*v.s + v.tau*v.kappa;
        
        %v.lam = full(sqrt(v.dxc'*(v.F{3}*v.dxc) ) );
        % same as: (cheaper to compute)
        v.lam = sqrt(v.dxc'*(rs.r5-v.dsc)/v.mu);
        
        % for printing:
        R.trace.cent.a(v.i)         = v.ac;
        R.trace.cent.lam(v.i)       = v.lam;
        R.trace.cent.tau(v.i)       = v.tau;
        R.trace.cent.kappa(v.i)     = v.kappa;
        R.trace.cent.cobj(v.i)      = v.centobjval;
        R.trace.cent.Pres(v.i)      = v.rPrel;
        R.trace.cent.Dres(v.i)      = v.rDrel;
        R.trace.cent.Gres(v.i)      = v.rGrel;
        R.trace.cent.Ares(v.i)      = v.rArel;
        R.trace.cent.by(v.i)        = pars.b'*v.y;
        R.trace.cent.cx(v.i)        = pars.c'*v.x;
        
        if v.lam < pars.beta
            break;
        end   
        
    end
    if i == pars.innermaxit
        v.reachedinnermaxit = v.reachedinnermaxit + 1;
    end
    if v.reachedinnermaxit > 2
        %fprintf('Reached innermaxit 3 times. Aborting \n');
        R.stop = 1;
        R.status = 'innermaxit3';
    end
    
end

R.trace.centstps(v.k+1)    = fullcentsteps;
R.trace.nbfgsstepst(v.k+1) = tbfgssteps;







