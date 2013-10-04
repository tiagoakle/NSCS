%this version is the one im making again
function R = coneopt(varargin)

[v,K,pars,R]  = init(varargin{:});

addpath '../C/interface'  %Add the diretory with the interface m files
addpath '../C/test'
load_c_library()
key = true; %Debug

while true
    
    if R.stop, break; end                % if R.stop, break loop
    
    v.F   = BarrFuncP(v.x,K,[1,1,1]);    % evaluate barrier
    
    %--------------------- Solve for the tangent direction ----------------------
    %this function implements the solve in the C function
    %[d,CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...
    % -v.rP,-v.rD,-v.rG,-v.kappa*v.tau,-v.s,pars);

    [m,n] = size(pars.A);

   % K5= [[sparse(m,m), pars.A             ,-pars.b    ,sparse(m,n),sparse(m,1)];...
   %      [-pars.A'   , sparse(n,n)        ,pars.c     ,-speye(n,n),sparse(n,1)];...
   %      [pars.b'    , -pars.c'           ,sparse(1,1),sparse(1,n),-1         ];...
   %      [sparse(1,m), sparse(1,n)        ,v.kappa    ,sparse(1,n),v.tau      ];...
   %      [sparse(n,m), sparse(v.mu*v.F{3}),sparse(n,1),speye(n,n) ,sparse(n,1)]];
   % 
   % rhs     = [-v.rP;-v.rD;-v.rG;-v.kappa*v.tau;-v.s];
   % d       = K5\rhs;
   % v.dy    = d(1:m);
   % v.dx    = d(m+1:m+n);
   % v.dtau  = d(m+n+1);
   % v.ds    = d(m+n+2:m+2*n+1);
   % v.dkappa= d(m+2*n+2);
    
    [d,CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,-v.rP,-v.rD,-v.rG,-v.kappa*v.tau,-v.s,pars);
    
   % fprintf('Differences '); 
   % vars = {'dx','dtau','dy','ds','dkappa'};
   % for j=1:5
   %   diff = eval(sprintf('norm( v.%s -d_2{%i})',vars{j},j));
   %   fprintf(' %s %g, ',vars{j},diff);
   % end
   % fprintf('\n');
   
   
   % counting:
    R.dat.nkktsolves = R.dat.nkktsolves  + 1;

   % if d{6} > 0
   %     R.status = 'rounding (atd)';
   %     R.stop   = true;
   % end

    v.dx     = d{1};
    v.dtau   = d{2};
    v.dy     = d{3};
    v.ds     = d{4};
    v.dkappa = d{5};

    v.dtauaff   = v.dtau;
    v.dkappaaff = v.dkappa;

    if R.stop                            % stop if rounding
        R = checkstopcrit(v,pars,R);     % blocks further progress
        break;                              
    end             
    %-------------------- End of tangent solve ----------------------------------
        
    %-------------------- Start of backtracking linesearch ----------------------
    %Call the C linesearch
    c_a = 1.0;
%    if v.k == 1
    [c_a,nbisections] = line_search_c(v,K,pars);
    v.a = c_a;
    xa = v.x + c_a*v.dx;
    sa = v.s + c_a*v.ds;
    taua = v.tau+c_a*v.dtau;
    kappaa = v.kappa+c_a*v.kappa;

%    v.a, nbisections
%    end

%    % set intial step length 
%    a0 = 1.0;
%    a  = a0;
%    if v.dkappa < 0
%        kapmax = -v.kappa/v.dkappa;
%        a = min(a,kapmax);
%    end
%    if v.dtau < 0
%        taumax = -v.tau/v.dtau;
%        a = min(a,taumax);
%    end
%    % if either kap or tau is blocking, multiply by eta
%    % so we do not hit boundary:
%    if a < a0
%        a = pars.eta*a;
%    end
%
%    % couter for number of bisections
%    nsect = 0;
%    
%    %Main linesearch loop
%    for j = 1:pars.lsmaxit
%        
%        % point to try:
%        xa     = v.x     + a * v.dx;
%        sa     = v.s     + a * v.ds;
%        taua   = v.tau   + a * v.dtau;
%        kappaa = v.kappa + a * v.dkappa;
%        
%        % new duality gap:
%        dga    = xa'*sa + taua*kappaa;
%        mua    = dga / (K.nu + 1);
%        
%        % evaluate barriers at new point:
%        % check only feasibility, so want = [-1,-1,-1]:
%        % Evaluate the primal barrier for f,g,H
%        FP = BarrFuncP(xa,K,[1,1,1]); 
%
%        %Call the c Hessian evaluation
%        primal_feasible = eval_primal_feas_c(K,xa);
%
%        
%        if(~isempty(FP{3})) 
%            [H] = eval_hessian_c(K,xa);
%            g   = eval_grad_c(K,xa);
%            HD  = full(max(max(abs(H-FP{3}))));
%            GD  = max(abs(g-FP{2}));
%            if HD > 1e-5
%                fprintf(' LARGE Hessian difference %g ', HD );
%            end
%            if GD > 1e-5
%                fprintf('LARGE Gradient difference %g ', GD );
%            end
%        end
%        % Check the dual feasibility
%        FD = BarrFuncD(sa,K,[1,-1,-1]);
%        dual_feasible = eval_dual_feas_c(K,sa);
%
%        if(FP{4}~=primal_feasible)
%            fprintf('Inconsistent feasiblity criteria for the primal (C:%i,M:%i)\n',primal_feasible,FP{4});
%        end
%        if(dual_feasible~=FD{4})
%            fprintf('Inconsistent feasiblity criteria for the dual (C:%i,M:%i)\n',dual_feasible,FD{4});
%        end
%        
%
%        dosect = false; %True if we must backtrack
%        %If either the primal is infeasible or 
%        % the dual is infeasible backtrack
%        if ~dual_feasible < 0 
%            dosect  = true;
%            R.block = 'pf';
%        elseif ~primal_feasible < 0 
%            dosect  = true;
%            R.block = 'df'; 
%        else %If the iterate is pirmal and dual feasible evaluate the centrality
%            
%           % psi       = sa + mua*FP{2};
%           % centmeas5 = sqrt(psi'*(FP{3}\psi)); %XXX: Linear solve
%           % centmeas = centmeas5;
%            
%            [centmeas] = eval_centmeas_c(K,xa,sa,mua);
%           
%         %   %XXX debug
%         %   m_hpsi = FP{3}\psi;
%         %   fprintf('Ndiff of Hinvpsi: %g\n', norm(hpsi-m_hpsi)); 
%         %   fprintf('Ndiff of psi: %g\n', norm(psi-c_psi)); 
%         %   fprintf('Ndiff of H %g\n',full(max(max(abs(H-FP{3})))));
%         %   fprintf('Ndiff of g %g\n',full(max(abs(gr-FP{2}))));
%         %   fprintf('mua_mat %f\n',mua);
%         %   fprintf('Centmeas: C: %f M: %f, lim: %f\n',c_centmeas, centmeas, mua*pars.theta);
%
%            if centmeas > mua*pars.theta
%                dosect  = true;
%                R.block = 'ce';
%            end
%        end
%        
%        if dosect
%            a     = a*pars.lscaff; 
%            nsect = nsect + 1;
%        else
%            break;
%        end
%            
%    end %end of main linesearch loop
%
%    v.a = a;
%    fprintf('C linesearch a, M linesearch a: %g, %g \n',c_a,v.a);
%    
%    % Check if the linesearch did not find 
%    % a feasible point in the maximum number of iterations
%    if j == pars.lsmaxit
%        xa = v.x + a * v.dx;
%        FP = BarrFuncP(xa,K,[1,1,-1]);
%        if FP{4} < 0
%
%            error(['linesearch: failed to find feasible point.',...
%                ' pars.lscaff too close to 1 ???']);
%        end
%    end

    %Take the step 

    % store the previous step
%    v.xprev     = v.x;
%    v.gprev     = v.F{2};
     v.tauprev   = v.tau;
     v.kappaprev = v.kappa;
    
    % take step:
    v.x     = xa;
    v.tau   = taua;
    v.y     = v.y     + v.a * v.dy;
    v.s     = sa;
    v.kappa = kappaa;
    
    % update other quantities:
    v.dgap  = v.x'*v.s + v.tau*v.kappa;
    v.mu    = v.dgap / (K.nu + 1);
    
    % "feas" measure, see Sturm:
    v.feas  = v.dtauaff/v.tauprev - v.dkappaaff/v.kappaprev;

   % %Update the residuals
   % bty  = pars.b'*v.y;
   % ctx  = pars.c'*v.x;    
   % v.rA = abs( ctx - bty )/( v.tau + abs(bty) );
   % 
   % v.rP = (pars.A*v.x - pars.b*v.tau);
   % v.rD = (-pars.A'*v.y - v.s + pars.c*v.tau);
   % v.rG = (-ctx + bty - v.kappa);
   % 
   % v.rPrel = norm( v.rP, 'inf')/pars.relstopP;
   % v.rDrel = norm( v.rD, 'inf')/pars.relstopD;
   % v.rGrel = norm( v.rG, 'inf')/pars.relstopG;
   % v.rArel = v.rA; 
    
    %Evaluate the c residuals function
    [p_res,d_res,g_res,n_p_res,n_d_res,n_g_res,rel_gap] = build_residuals_c(K,v.y,v.x,v.s,v.tau,v.kappa,pars.A,pars.b,pars.c,pars.relstopP,pars.relstopD,pars.relstopG);

    %Update the residuals   
    v.rA = rel_gap;
    
    v.rP = -p_res;
    v.rD = -d_res;
    v.rG = -g_res;
    
    v.rPrel = n_p_res;
    v.rDrel = n_d_res;
    v.rGrel = n_g_res;
    v.rArel = v.rA; 
   
    
    %------------------------ End of backtracking linesearch  ------------

    R       = checkstopcrit(v,pars,R);   % check stopping crits
    
    if R.stop
        fprintf('Stop command issued\n');    
    end

    %------------------------- Centering process -----------------------

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
            keyboard
            rs.r5 = - v.s-v.mu*v.F{2}; %        - (-pars.A'*v.y + v.mu*v.F{2} + v.tau*pars.c);
         
            % compute the centering measure at the present point
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

%            K5= [[sparse(m,m), pars.A             ,-pars.b    ,sparse(m,n),sparse(m,1)];...
%                 [-pars.A'   , sparse(n,n)        ,pars.c     ,-speye(n,n),sparse(n,1)];...
%                 [pars.b'    , -pars.c'           ,sparse(1,1),sparse(1,n),-1         ];...
%                 [sparse(1,m), sparse(1,n)        ,v.kappa    ,sparse(1,n),v.tau      ];...
%                 [sparse(n,m), sparse(v.mu*v.F{3}),sparse(n,1),speye(n,n) ,sparse(n,1)]];
%            
%            rhs     = [rs.r1;rs.r2;rs.r3;rs.r4;rs.r5];
%            d       = K5\rhs;
%            v.dyc    = d(1:m);
%            v.dxc    = d(m+1:m+n);
%            v.dtauc  = d(m+n+1);
%            v.dsc    = d(m+n+2:m+2*n+1);
%            v.dkappac= d(m+2*n+2);

            [d,CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,rs.r1,rs.r2,rs.r3,rs.r4,rs.r5);
            v.dxc    = d{1};
            v.dtauc  = d{2};
            v.dyc    = d{3};
            v.dsc    = d{4};
            v.dkappac= d{5};
 
            % counting:
            R.dat.nkktsolves = R.dat.nkktsolves  + 1;

            % stop if rounding is blocking progress:
            if R.stop, R = checkstopcrit(v,pars,R); break; end
            
            % line search:
            %v = linesearchcent(v,K,pars);
            [a_cent,nbacktracks,objvalb] = line_search_cent_c(v,K,pars);
            v.ac = a_cent;
            v.centobjval = objvalb;

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

    %-------------------------End of centering process--------------------

    v.k     = v.k + 1;                   % increment counter
    R       = gettrace(v,pars,R);        % get current data
    
    printfunc('iter',v,K,pars,R);        % print iteration info
    
end

[v,K,pars,R]  = deinit(v,K,pars,R);
printfunc('final',v,K,pars,R);

unload_c_library()
