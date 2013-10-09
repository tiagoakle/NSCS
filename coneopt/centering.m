function [v,R] = centering(v,K,pars,R)

v.i     = 0;
v.lam   = 99;
v.jt    = 0;


fullcentsteps = 0;
tbfgssteps    = 0;
stopcentering = false;

if ~R.stop
    
    v.rD0  = -pars.A'*v.y - v.s + v.tau*pars.c;
    
    for i = 1:pars.innermaxit
        v.i = i;
        
        % evalute barrier at new x:
        v.F = BarrFuncP(v.x,K,[1,1,1]);
        
        % right-hand sides:
        rs.r1 = v.rP   - (pars.A*v.x - pars.b*v.tau);
        rs.r2 = v.rD0;
        rs.r3 = v.rG   - (-pars.c'*v.x + pars.b'*v.y - v.kappa);
        rs.r4 = v.mu   - (v.tau*v.kappa);
        rs.r5 =        - (-pars.A'*v.y + v.mu*v.F{2} + v.tau*pars.c);
        
        
        % compute first centering measure:
        r2cent     = rs.r2 + rs.r5;
        tmp        = r2cent/v.mu;
        [L,pp]     = chol(v.F{3});
        if pp > 0
            v.centres  = 99;
        else
            tmp        = L'\tmp;
            v.centres  = sqrt(tmp'*tmp);
        end
        
        % Take pars.cnbfgsstps BFGS steps:
        v.Q1(:)   = 0;
        v.Q2(:)   = 0;
        v.MQ2(:)  = 0;
        v.sc1(:)  = 0;
        
        bfgssteps = 0;
        
        for j = 1:pars.cnbfgsstps
            bfgssteps  = bfgssteps + 1;
            tbfgssteps = tbfgssteps + 1;
            v.j = j;
            
            % compute ss and yy:
            ss = v.x - v.xprev;
            yy = v.mu*(v.F{2} - v.gprev);
            % notice: doing BFGS on mu*F, not F
            
            % check the curvature condition:
            if ss'*yy <= 0
                break;
            end
            
            % By = B*yy
            By = (v.CF.L)\((v.CF.L)'\yy) + v.Q1*(v.sc1.*(v.Q1'*yy));
            
            % BFGS update:
            upd = bfgsupd(ss,yy,By);
            
            % solve for new direction:
            [d,v] = slvbfgsdir_coneopt(v.CF,v,rs,upd,pars);
            
            % extract sol:
            v.dxc = d{1}; 
            v.dtauc = d{2}; 
            v.dyc = d{3};
            v.dsc = d{4}; 
            v.dkappac = d{5};
            
            % line search
            v = linesearchcent(v,K,pars);
            if v.ac < pars.stpszmin
                break;
            end
            
            % for next bfgs:
            v.xprev = v.x;
            v.gprev = v.F{2};
            
            % take step:
            v.x     = v.x     + v.ac * v.dxc;
            v.tau   = v.tau   + v.ac * v.dtauc;
            v.y     = v.y     + v.ac * v.dyc;
            v.kappa = v.kappa + v.ac * v.dkappac;
            v.s     = v.s     + v.ac * v.dsc;
            
            % new barrier:
            v.F = BarrFuncP(v.x,K,[1,1,1]);
            
            % new right-hand sides:
            rs.r1 = v.rP   - (pars.A*v.x - pars.b*v.tau);
            rs.r2 = v.rD0;
            rs.r3 = v.rG   - (-pars.c'*v.x + pars.b'*v.y - v.kappa);
            rs.r4 = v.mu   - (v.tau*v.kappa);
            rs.r5 = v.r5cent;
            
            % new centering measure:
            v.centres  = v.centobjval;
            
            % for printing:
            R.trace.bfgsc.a{v.i}(v.j)     = v.ac;
            R.trace.bfgsc.lam{v.i}(v.j)   = v.lam;
            R.trace.bfgsc.lamh{v.i}(v.j)  = v.lam;
            R.trace.bfgsc.tau{v.i}(v.j)   = v.tau;
            R.trace.bfgsc.kappa{v.i}(v.j) = v.kappa;
            R.trace.bfgsc.cobj{v.i}(v.j)  = v.centobjval;
            R.trace.bfgsc.Pres{v.i}(v.j)  = v.rPrel;
            R.trace.bfgsc.Dres{v.i}(v.j)  = v.rDrel;
            R.trace.bfgsc.Gres{v.i}(v.j)  = v.rGrel;
            R.trace.bfgsc.Ares{v.i}(v.j)  = v.rArel;
            
            if pars.bfgsstop && v.centres < pars.beta
                stopcentering = true;
                break;
            end
            
        end % END BFGS
        
        R.trace.cent.bfgsstps(v.i)  = bfgssteps;
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







