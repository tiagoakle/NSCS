%Minimal implementation of a mehrtra type predictor corrector for the non symmetric strategy
function [x,y,s,info] = non_symmetric_long_step(A,b,c,strategy)

    warning('off','MATLAB:nearlySingularMatrix');

    fprintf('Minimal Nonsymmetric mehrotra type solver \n');
    %-----------------------------------------------------------------
    % Minimal Mehrotra LP predictor corrector
    %------------------------------------------------------------------
    [m,n] = size(A); 
    nu    = n;
    fprintf('Problem size %i constraints %i variables \n',m,n);

    %-----------------------------------
    %Initialization strategy from CVXOPT
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[zeros(n,1);b];
    x   = sol(1:n);
    
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[c;zeros(m,1)];
    y   = sol(n+1:n+m);
    s   = sol(1:n);
    
    %%XXX:RANDOM INIT
    %x   = ones(n,1);
    %s   = ones(n,1);
    %y   = zeros(m,1);
    
    tau   = 1;
    kappa =1;
    clear sol
    
    if(min(x)<0) %if x is not feasible shift it into feasibility
        a = min(x);
        x = (1-a)*ones(n,1)+x;
    end
    
    if(min(s)<0) %if s is not feasible shift it into feasibility
        a = min(s);
        s = (1-a)*ones(n,1)+s;
    end
    
    %------Calculate the initial mu
    mu = (s'*x+tau*kappa)/(nu+1);
    mu0 = mu;
    
    %--------------------
    %Evaluate the residuals
    rp = tau*b-A*x;
    rd = -tau*c+A'*y+s;
    rg = kappa-b'*y+c'*x;
    
    %Evaluate the residual norms
    nrp = norm(rp);
    nrd = norm(rd);
    nrg = norm(rg);
    
    % Save the original norms
    nrp0 = nrp;
    nrd0 = nrd;
    nrg0 = nrg;    
    gap  = c'*x-b'*y;

    %Algorithm parameters
    max_arc_backtrack_iter = 300;
    arc_search_backtrack_scale = 0.8; 
    eta_scaling            = 0.8;
    max_iter = 100;

    fprintf('%2i a       %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',0,nan,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,nan);
    exit_reason = '';

    a = 1;
    %--------------------------------------
    %Main iteration
    for iter=1:max_iter
         
        %Evaluate the hessian and scaled variables 
        %Solve the aproximate affine scaling direction
        H = diag(sparse(mu*(1./x).^2));
        [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-s,-kappa,[]);

        ratios = [-x./dx;-s./ds;-tau/dtau;-kappa/dkappa;1];
        a_max  = min(ratios(find(ratios>0)));
    
        %Now calculate the centering direction
        %-------------------------------------
        sigma = (1-a_max)^3;
        
        %These variables hold the second order info
        so = 0;
        kt_so = 0;
        if(~strcmp(strategy,'linear'))
            temp = s+ds;
            so   =  (a)*(0.5/mu*x.*(temp.^2)-temp);
            kt_so= -(a)*dtau*dkappa/tau;
            clear temp
        end

        if(strcmp(strategy,'arc_search'))
            %Solve for the second derivative of the path
            [dy_s,dx_s,dtau_s,ds_s,dkappa_s,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,...
                                                                                                zeros(m,1),...
                                                                                                zeros(n,1),...
                                                                                                0,...
                                                                                                so,...
                                                                                                kt_so,...
                                                                                                slv_aug);
            so = 0;
            kt_so = 0; 
        end
        %Solve the affine scaling direction
        [dy_c,dx_c,dtau_c,ds_c,dkappa_c,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,...
                                                                                                (1-sigma)*rp,...
                                                                                                (1-sigma)*rd,...
                                                                                                (1-sigma)*rg,...
                                                                                                -s    +sigma*mu*1./x  + (1-sigma)*so,...
                                                                                                -kappa+sigma*mu*1/tau + (1-sigma)*kt_so,...
                                                                                                slv_aug);
    
        if(~strcmp(strategy,'arc_search')) %Mehrotra type step for secord and no secord
            %Find the largest step to the boundary
            ratios = [-x./dx_c;-s./ds_c;-tau/dtau_c;-kappa/dkappa_c;1];
            a_max  = min(ratios(find(ratios>0)));
            a       = a_max*eta_scaling; 
            
            x_old = x;
            y_old = y;
            s_old = s;
            t_old = tau;
            k_old = kappa;
    
            y         = y+a*dy_c;
            x         = x+a*dx_c;
            s         = s+a*ds_c;
            tau       = tau+a*dtau_c;
            kappa     = kappa+a*dkappa_c;
        else %Do the arc search
            %Backtrack into feasibility
            a = 1;
            x_old = x;
            y_old = y;
            s_old = s;
            t_old = tau;
            k_old = kappa;
            
            for(j=1:max_arc_backtrack_iter)
                y         = y_old+a*dy_c+a^2*dy_s;
                x         = x_old+a*dx_c+a^2*dx_s;
                s         = s_old+a*ds_c+a^2*ds_s;
                tau       = t_old+a*dtau_c+a^2*dtau_s;
                kappa     = k_old+a*dkappa_c+a^2*dkappa_s;
                if((min(x)>0)&&(min(s)>0)&&(min(tau,kappa)>0))
                    break;
                end
                a = a*arc_search_backtrack_scale;
            end
            if(j==max_arc_backtrack_iter)
                fprintf('Unable to find feasible arc point\n');
                exit_reason = 'Backtrack fail';
                break;
            end
            if(j>1)
                %one more scaling to move away from the boundary
                a = a*eta_scaling;
            end

            y         = y_old+a*dy_c+a^2*dy_s;
            x         = x_old+a*dx_c+a^2*dx_s;
            s         = s_old+a*ds_c+a^2*ds_s;
            tau       = t_old+a*dtau_c+a^2*dtau_s;
            kappa     = k_old+a*dkappa_c+a^2*dkappa_s;

        end

        %Update mu
        mu        = x'*s+tau*kappa;
        mu        = mu/(nu+1);
        
        %Feasiblity check
        x_slack = min(x);
        s_slack = min(s);
        if(x_slack<0||s_slack<0||tau<0||kappa<0)
            fprintf('Infeasible centering step, a: %g \n',a);
            exit_reason = 'Infeas centering';
            a = 1;
            min(x_old+a*dx_c+a^2*dx_s)
            break;
        end
    
        %--------------------
        %Evaluate the residuals
        rp = tau*b-A*x;
        rd = -tau*c+A'*y+s;
        rg = kappa-b'*y+c'*x;
    
        nrp   = norm(rp);
        nrd   = norm(rd);
        nrg   = norm(rg);
        
        gap   = (mu*(nu+1)-tau*kappa)/tau;
        
        fprintf('%2i a %3.3e s %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',iter,a,sigma,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,residual_norm_c);
        if(nrp/(nrp0) < 1.e-8 && nrd/(nrd0) < 1.e-8 && mu/mu0 < 1.e-8)
            exit_reason = 'finished';
            break;
        end
    end %end main loop

    x = x/tau;
    s = s/tau;
    y = y/tau;
    info = struct;
    info.iter = iter;
    if(~strcmp(exit_reason,'finished'))
        info.iter = -1;
    end
end

