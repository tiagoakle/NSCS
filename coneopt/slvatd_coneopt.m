function [v,R] = slvatd_coneopt(v,K,pars,R)

% solve for
% approximate tangent direction

r4 = -v.kappa*v.tau;
r5 = -v.s;

% another option that could be investigated:
%r5 = v.mu * v.F{2};

[d,CF] = slvhomkkt(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
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



% TESTING
if pars.addcent
beta1 = 0.9;
[am,R] = linesearch(v,K,pars,R);
gam    = (1-am)*min((1-am)^2,beta1);
    
% solve for final search dir:
r1  = -v.rP          + gam*v.rP;
r2  = -v.rD          + gam*v.rD;
r3  = -v.rG          + gam*v.rG;
r4  = -v.tau*v.kappa + gam*v.mu;
r5  = -v.s           - gam*v.mu*v.F{2};

[d,CF] = slvhomkkt(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
    r1,r2,r3,r4,r5,pars);
    
v.dx     = d{1};
v.dtau   = d{2};
v.dy     = d{3};
v.ds     = d{4};
v.dkappa = d{5};
end




v.CF = CF;

% SECORD1 (real Hessian at xi)
if pars.secord
    
    % determine h
    [h,R] = linesearch(v,K,pars,R);
    
    if h < pars.hlim
        
        p     = pars.RK2param;
        a21   = p;
        b2    = 1/(2*p);
        b1    = 1-b2;
        ha    = a21*h;
        
        % compute xi:
        xix     = v.x     + ha * v.dx;
        xis     = v.s     + ha * v.ds;
        xiy     = v.y     + ha * v.dy;
        xitau   = v.tau   + ha * v.dtau;
        xikappa = v.kappa + ha * v.dkappa;
        
        % update at xi:
        dgxi    = xix'*xis + xitau*xikappa;
        muxi    = dgxi / (K.nu + 1);
        Fxi     = BarrFuncP(xix,K,[1,1,1]);
        rPxi    = (pars.A*xix - pars.b*xitau);
        rDxi    = (-pars.A'*xiy - xis + pars.c*xitau);
        rGxi    = (-pars.c'*xix + pars.b'*xiy - xikappa);
        
        r4 = -xikappa*xitau;
        if pars.inclmu
            r4 = r4 + muxi;
        end
        
        % solve for new direction:
        d2 = slvhomkkt(Fxi{3},muxi,pars.A,pars.b,pars.c,xitau,xikappa,...
            -rPxi,-rDxi,-rGxi,r4,-xis,pars);
        
        % counting:
        R.dat.nkktsolves = R.dat.nkktsolves  + 1;
        
        if d2{6} > 0
            R.status = 'rounding (atd)';
            R.stop   = true;
            return;
        end
        
        % new total search dir:
        for j = 1:5
            d{j} = 2*h*(b1*d{j} + b2*d2{j});
        end
        
        v.dx     = d{1};
        v.dtau   = d{2};
        v.dy     = d{3};
        v.ds     = d{4};
        v.dkappa = d{5};
        
    end
    
% SECORD2 (BFGS at xi)
elseif pars.secord2
    
    
    % determine h
    [h,R] = linesearch(v,K,pars,R);
    
    if h < pars.hlim
        
        p     = pars.RK2param;
        a21   = p;
        b2    = 1/(2*p);
        b1    = 1-b2;
        ha    = a21*h;
        
        % compute xi:
        xix     = v.x     + ha * v.dx;
        xis     = v.s     + ha * v.ds;
        xiy     = v.y     + ha * v.dy;
        xitau   = v.tau   + ha * v.dtau;
        xikappa = v.kappa + ha * v.dkappa;
        
        % update at xi:
        dgxi    = xix'*xis + xitau*xikappa;
        muxi    = dgxi / (K.nu + 1);
        Fxi     = BarrFuncP(xix,K,[1,1,1]);
        rPxi    = (pars.A*xix - pars.b*xitau);
        rDxi    = (-pars.A'*xiy - xis + pars.c*xitau);
        rGxi    = (-pars.c'*xix + pars.b'*xiy - xikappa);
        
        % right hand sides:
        r1 = -rPxi;
        r2 = -rDxi;
        r3 = -rGxi;
        r4 = -xikappa*xitau;
        r5 = -xis;
        
        % faster BFGS: -------------------
        Q1  = zeros(pars.n,2);
        Q2  = zeros(pars.m,2);
        MQ2 = zeros(pars.m,2);
        sc1 = zeros(2,1);
        L   = CF.L;
        N   = CF.N;
        
        % compute ss and yy:
        ss = xix   - v.x;
        yy = v.mu*(Fxi{2} - v.F{2});
        % notice: doing BFGS on mu*F, not F
        
        % check the curvature condition:
        if ss'*yy <= 0
            error('atd: bfgs direction bad')
        end
        
        % By = B*yy
        By = L\(L'\yy);
        
        % BFGS update:
        upd = bfgsupd(ss,yy,By);
        
        % update the BFGS data holders:
        idx1             = 1;
        idx2             = 2;
        Q1(:,idx1:idx2)  = [upd.v,upd.w];
        sc1(idx1:idx2)   = [upd.sgnv;upd.sgnw];
        Q2(:,idx1:idx2)  = pars.A*Q1(:,idx1:idx2);
        if pars.permuteM
            tmp              = pars.Sp'*Q2(:,idx1:idx2);
            tmp              = N\(N'\tmp);
            MQ2(:,idx1:idx2) = pars.Sp*tmp;
        else
            MQ2(:,idx1:idx2) = N\(N'\(Q2(:,idx1:idx2)));
        end
        Q3               = Q2(:,1:idx2);
        MQ3              = MQ2(:,1:idx2);
        sc3              = sc1(1:idx2);
        MS               = (diag(sc3) + Q3'*MQ3);
        
        % solve sequence:=======
        sc1r  = repmat(sc1,1,2);
        r6    = r2 + r5;
        tmp   = [-r6,pars.c];
        % tmp = B*tmp
        tmp   = L\(L'\tmp) + Q1*(sc1r.*(Q1'*tmp));
        tmp   = [r1,pars.b] + pars.A*tmp;
        % this is tmp = M*tmp, where M ~= (A*H^{-1}*A')^{-1}
        tmp   = tmp - Q3*( MS\(MQ3'*tmp) );
        if pars.permuteM
            tmp = pars.Sp'*tmp;
        end
        tmp   = N\(N'\tmp);
        if pars.permuteM
            tmp = pars.Sp*tmp;
        end
        % --------------------
        r9    = tmp(:,1);
        r10   = tmp(:,2);
        tmp   = [r6,-pars.c] + pars.A'*tmp;
        % tmp = B*tmp
        tmp   = L\(L'\tmp) + Q1*(sc1r.*(Q1'*tmp));
        r11   = tmp(:,1);
        r12   = tmp(:,2);
        % for ease of notation (and copying later):
        c     = pars.c;
        b     = pars.b;
        A     = pars.A;
        tau   = xitau;
        kappa = xikappa;
        % ------------------
        t1    = r4 + tau*(c'*r11 - b'*r9 + r3);
        t2    = tau*(b'*r10 - c'*r12) + kappa;
        d2{2}  = t1/t2;
        d2{5}  = (r4 - kappa*d2{2})/tau;
        d2{3}  = r9  + r10*d2{2};
        d2{1}  = r11 + r12*d2{2};
        d2{4}  = c*d2{2} - r2 - A'*d2{3};
        % == solve end ============
        
        % new total search dir:
        for j = 1:5
            d{j} = 2*h*(b1*d{j} + b2*d2{j});
        end
              
        v.dx     = d{1};
        v.dtau   = d{2};
        v.dy     = d{3};
        v.ds     = d{4};
        v.dkappa = d{5};
        
    end
    
% test adding some centering to the right hand side:
elseif pars.secord3
    
    beta1  = 0.1;
    [am,R] = linesearch(v,K,pars,R);
    gam    = (1-am)^2*min((1-am),beta1);
    
    % solve for final search dir:
    r1  = -v.rP          + gam*v.rP;
    r2  = -v.rD          + gam*v.rD;
    r3  = -v.rG          + gam*v.rG;
    r4  = -v.tau*v.kappa + gam*v.mu;
    r5  = -v.s           - gam*v.mu*v.F{2};
    
    [d,CF] = slvhomkkt(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
        r1,r2,r3,r4,r5,pars);
    
    v.dx     = d{1};
    v.dtau   = d{2};
    v.dy     = d{3};
    v.ds     = d{4};
    v.dkappa = d{5};
    
    
    
elseif pars.secord4
    
     % determine h
    [h,R] = linesearch(v,K,pars,R);
    
    if h < pars.hlim
        
        p     = pars.RK2param;
        a21   = p;
        b2    = 1/(2*p);
        b1    = 1-b2;
        ha    = a21*h;
        
        % compute xi:
        xix     = v.x     + ha * v.dx;
        xis     = v.s     + ha * v.ds;
        xiy     = v.y     + ha * v.dy;
        xitau   = v.tau   + ha * v.dtau;
        xikappa = v.kappa + ha * v.dkappa;
        
        % update at xi:
        dgxi    = xix'*xis + xitau*xikappa;
        muxi    = dgxi / (K.nu + 1);
        Fxi     = BarrFuncP(xix,K,[1,1,1]);
        
        %rPxi    = (pars.A*xix - pars.b*xitau);
        %rDxi    = (-pars.A'*xiy - xis + pars.c*xitau);
        %rGxi    = (-pars.c'*xix + pars.b'*xiy - xikappa);       
        
        % compute ss and yy:
        ss = xix   - v.x;
        yy = v.mu*(Fxi{2} - v.F{2});
        % notice: doing BFGS on mu*F, not F
        
        % check the curvature condition:
        if ss'*yy <= 0
            %do something here
        end
        
        % By = B*yy
        By = CF.L\(CF.L'\yy);
        
        % BFGS update:
        upd = bfgsupd(ss,yy,By);
        
        v.Q1(:)  = 0;
        v.Q2(:)  = 0;
        v.MQ2(:) = 0;
        v.sc1(:) = 0;
        jp = v.j; v.j = 1;
        rs.r1  = zeros(pars.m,1);
        rs.r2  = zeros(pars.n,1);
        rs.r3  = zeros(1,1);
        rs.r4  = -xikappa*xitau + muxi;
        rs.r5  = -xis - muxi*Fxi{2};
        
        [d2,v] = slvbfgsdir_coneopt(CF,v,rs,upd,pars);
        v.j    = jp;
        
        % new total search dir:
        for j = 1:5
            d{j} = 4*h*(b1*d{j} + b2*d2{j});
        end
        
        v.dx     = d{1};
        v.dtau   = d{2};
        v.dy     = d{3};
        v.ds     = d{4};
        v.dkappa = d{5};
        
    end
    
    
end






