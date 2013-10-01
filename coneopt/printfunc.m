function printfunc(s,v,K,pars,R)

% echo:
% 0  = nothing at all
% 1  = limited output
% 2  = narrow, only outer iter
% 3  = wide1, only outer iter
% 4  = wide2, only outer iter
% 12 = as 2, but with inner iter
% 13 = as 3, but with inner iter
% 14 = as 4, but with inner iter
% -x = as x but all printed post-run, nothing during run.


if pars.echo > 0
    
    switch s
        
        % OUTER ITERATION -----------------------------------
        case 'iter'
            
            it          = v.k-1;
            ncentstepe  = R.trace.ncentstepse(v.k);
            nbfgsstepe  = R.trace.nbfgsstepst(v.k);
            fcentsteps  = R.trace.centstps(v.k);
            a           = R.trace.a(v.k);
            tau         = R.trace.tau(v.k);
            kap         = R.trace.kappa(v.k);
            mu          = R.trace.mu(v.k);
            Ares        = R.trace.rArel(v.k);
            Gres        = R.trace.rGrel(v.k);
            Pres        = R.trace.rPrel(v.k);
            Dres        = R.trace.rDrel(v.k);
            feas        = R.trace.feas(v.k);
            PDres       = max(Pres,Dres);
            AGres       = max(Ares,Gres);
            cx          = pars.c'*v.x;
            by          = pars.b'*v.y;
            dgap        = v.dgap;
            
            
            if pars.echo > 10, printline(pars.echo,1), end;
            
            if rem(pars.echo,10) == 2
                fprintf(...
                    [' %-3.1i %-3.1i %-3.1i %9.2E %9.2E ',...
                    '%9.2E %9.2E %9.2E %9.2E',...
                    '\n'],...
                    it, fcentsteps, nbfgsstepe, ...
                    a, mu, kap/tau, feas, PDres, AGres );
                
            elseif rem(pars.echo,10) == 3
                fprintf(...
                    [' %-3.1i %-3.1i %9.2E %9.2E ',...
                    '%9.2E %9.2E %9.2E',...
                    ' %+9.7E %+9.7E \n'],...
                    it, fcentsteps, a, mu,...
                    kap/tau,PDres,AGres,cx/tau,by/tau);
                
            elseif rem(pars.echo,10) == 4
                fprintf(...
                    [' %-3.1i %-3.1i %9.2E %9.2E ',...
                    '%9.2E %9.2E ',...
                    '%9.2E %9.2E %9.2E ',...
                    '%9.2E \n'],...
                    it, fcentsteps, a, mu,...
                    tau,kap,Pres, Dres, Ares,dgap);
            end
            
            if pars.echo > 10, printline(pars.echo,1), end;
            % END OUTER ITER -----------------------------------------
            
            % INNER ITERS ----------------------------
            if abs(pars.echo) > 10
                for i = 1:ncentstepe
                    
                    bfgsstps = R.trace.cent.bfgsstps(i);
                    
                    % BFGS iters:
                    if abs(pars.echo) > 20
                        for j = 1:bfgsstps
                            
                            ac      = R.trace.bfgsc.a{i}(j);
                            lam     = R.trace.bfgsc.lam{i}(j);
                            lamh    = R.trace.bfgsc.lamh{i}(j);
                            tau     = R.trace.bfgsc.tau{i}(j);
                            kappa   = R.trace.bfgsc.kappa{i}(j);
                            cobj    = R.trace.bfgsc.cobj{i}(j);
                            Pres    = R.trace.bfgsc.Pres{i}(j);
                            Dres    = R.trace.bfgsc.Dres{i}(j);
                            Gres    = R.trace.bfgsc.Gres{i}(j);
                            Ares    = R.trace.bfgsc.Ares{i}(j);
                            mPDGres = max([Pres,Dres,Gres]);
                            
                            if rem(pars.echo-20,10) == 2
                                
                                fprintf(...
                                    ['         %-3.1i %9.2E %9.2E ',...
                                    '%9.2E %9.2E %9.2E \n'],...
                                    j, ac, lamh, kappa/tau, cobj, mPDGres);
                                
                            elseif rem(pars.echo-20,10) == 3
                                
                            elseif rem(pars.echo-20,10) == 4
                                
                                
                                
                                
                            end
                        end
                    end
                    
                    if abs(pars.echo) > 10
                        if fcentsteps > i-1
                            
                            ac       = R.trace.cent.a(i);
                            lam      = R.trace.cent.lam(i);
                            kap      = R.trace.cent.kappa(i);
                            tau      = R.trace.cent.tau(i);
                            cobj     = R.trace.cent.cobj(i);
                            rPrel    = R.trace.cent.Pres(i);
                            rDrel    = R.trace.cent.Dres(i);
                            rGrel    = R.trace.cent.Gres(i);
                            mPDGres  = max([rPrel,rDrel,rGrel]);
                            cx       = R.trace.cent.cx(i);
                            by       = R.trace.cent.by(i);
                            
                            
                            if rem(pars.echo-10,10) == 2
                                
                                fprintf(...
                                    ['     %-3.1i %-3.1i %9.2E %9.2E ',...
                                    '%9.2E %9.2E %9.2E \n'],...
                                    i, bfgsstps, ac, lam, ...
                                    kap/tau, cobj, mPDGres);
                                
                            elseif rem(pars.echo-10,10) == 3
                                
                                fprintf(...
                                    ['     %-3.1i %9.2E %9.2E ',...
                                    '%9.2E %9.2E %9.2E %+9.7E %+9.7E \n'],...
                                    i, ac, lam, kap/tau, cobj,...
                                    mPDGres,cx,by);
                                
                                
                            elseif rem(pars.echo-10,10) == 4
                                
                                fprintf(...
                                    ['     %-3.1i %9.2E %9.2E ',...
                                    '%9.2E %9.2E %9.2E %9.2E %9.2E %9.2E \n'],...
                                    i, ac, lam, tau, kap,...
                                    rPrel,...
                                    rDrel,...
                                    cobj,cx-by);
                                
                            end
                        end
                    end
                    
                    
                    
                    
                end
            end
            % END INNER ITER --------------------------------
            
            
            
            % INITIAL PRINTING (HEADER) ------------------
        case 'initial'
            printline(pars.echo,2);
            fprintf(' problem has %i variables and %i eq-constraints. \n',pars.n,pars.m);
            %fprintf(' Stopping: abs/rel dg-tol: %8.2e / %8.2e. \n',pars.eps*pars.relstopG,pars.eps);
            printline(pars.echo,2);
            if pars.echo < 2
                fprintf(' Running...\n');
            else
                if rem(pars.echo,10) == 2
                    fprintf([' oi  ii  bi   alpha     mu        kap/tau  ',...
                        ' feas      PD-res    GA-res   \n']);
                elseif rem(pars.echo,10) == 3
                    fprintf([' it  ii   alpha     mu        kap/tau  ',...
                        ' P&D-res   G&A-res   P-objval       D-objval \n']);
                elseif rem(pars.echo,10) == 4
                    fprintf([' it  ii   alpha     mu        tau       kappa  ',...
                        '   p-res     d-res     A-res     d-gap \n']);
                end
            end
            
            printline(pars.echo,1);
            if rem(pars.echo-10,10) == 2
                fprintf(['     ii  bi   a         lam       ',...
                    'kap/tau   c-obj     PDG-res \n']);
            elseif rem(pars.echo-10,10) == 3
                fprintf(['     ii   a         lam       ',...
                    'kap/tau   c-obj     PDG-res   ',...
                    'P-objval       D-objval \n']);
            elseif rem(pars.echo-10,10) == 4
                fprintf(['     ii   a         lam       tau       kappa',...
                    '     p-res     d-res     c-obj     d-gap \n']);
            end
            
            %printline(pars.echo,2);
            % END INITIAL (HEADER) ---------------------------------
            
            if pars.echo > 1
                printfunc('iter',v,K,pars,R);
            end
            
        case 'final'
        otherwise
            error('unknown type in printfunc');
    end
end


if abs(pars.echo) > 0
    if strcmp(s,'final')
        if pars.echo > 0
            printline(pars.echo,2);
            fprintf(' Terminated because: ');
            switch R.status
                case 'maxit'
                    str = 'Reached maxit';
                case 'optimal'
                    str = 'Optimal';
                case 'debug'
                    str = 'debug';
                case 'P-infeasible'
                    str = 'P-infeasible (D-unbounded)';
                case 'D-infeasible'
                    str = 'D-infeasible (P-unbounded)';
                case 'ill-posed'
                    str = 'Ill-posed problem';
                case 'rounding (atd)'
                    str = ['\n Rounding prevented further progress (atd).\n',...
                        ' rel-gap = ',num2str(v.mu,'%9.2E')];
                case 'rounding (cent)'
                    str = ['\n Rounding prevented further progress (cent).\n',...
                        ' rel-gap = ',num2str(v.mu,'%9.2E')];
                case 'innermaxit3'
                    str = 'innermaxit was reached 3 times. Aborting';
                otherwise
                    str = 'Something else happened';
            end
            fprintf([str,'\n']);
            printline(pars.echo,2);
            spcs = '          ';
            fprintf(' Used: %3.1e seconds \n',R.dat.tt);
            fprintf([spcs,fixedlengthint(v.k-1,4,'r'),...
                ' iterations \n']);
            fprintf([spcs,fixedlengthint(R.dat.ncentsteps,4,'r'),...
                ' centering steps \n']);
            fprintf([spcs,fixedlengthint(sum(R.trace.nbfgsstepst),4,'r'),...
                ' bfgs corrections \n']);
            fprintf([spcs,fixedlengthint(R.dat.nkktsolves,4,'r'),...
                ' kkt solves \n']);
            %fprintf([num2str(R.dat.tt,'%9.2E'),' seconds \n']);
            printline(pars.echo,2);
            % END FINAL PRINT -----------------------------
            
        else % pars.echo < 0, so print all in one go
            pars.echoorig = pars.echo;
            pars.echo     = -pars.echo;
            iters         = v.k;
            v.k           = 1;
            printfunc('initial',v,K,pars,R);
            for k = 2:iters
                v.k = k;
                printfunc('iter',v,K,pars,R);
                %if pars.echo > 10
                %    printfunc('inneriter',v,K,pars,R);
                %end
                pause(0.01);
            end
            printfunc('final',v,K,pars,R);
            
        end
    end
end






