%Sep 28 20x2
%Primal dual interior point Homogeneous self dual for non symmetric cones

%It saves the generated gradients and all
%iteration data.

%It solves CPs in standard form
% min c'x s.t. Ax=b and x in k

function [x,y,z,err_logs,iters_log] = prip_toy_nesterov(A,b,c,x0,t0,y0,z0,k0,pars)

addpath ../Solver
%Make sure we use this _minres
%Parameters
%A,b,c LP standard form parameters
%x0,y0,z0 initial point
%pars: parameter structure
%errs: error log structure

  %---------------------------------------
  %Define the parameters
  %---------------------------------------
  if(~exist('pars','var'))
    pars = struct;
  end
  
  
  % set defaults:
  defaults = {...
      'tol',1e-4;... %Tolerance to stop 
      'max_bkt',30;... %Backtracking iterations
      'beta_bkt',0.5;... %Backtracking constant
      'eta_bkt',0.8;...  %Fraction of decrease that must be achieved
      'max_iter',100;... %Iteration limit
      'max_compound',Inf;... %Maximum compound iterations
      'minres_maxit',Inf;...
      'rho',-1;...          %Rho -1 means that rho will be chosen as n+sqrt(n)
      'errFcn',{}; 
      'stopFcn',[];
  	  'save_iters',false;
      'tau',0.9;            %How close to the boundary to get at each step
   };
  
  for k = 1:size(defaults,1)
      if ~isfield(pars,defaults(k,1))
          pars.(defaults{k,1}) = defaults{k,2};
      end
  end

%Define the default return parameters
iters_log = [];
err_logs  = [];
[m,n] = size(A);

%Tolerances 
tol      = pars.tol;
max_iter = pars.max_iter;
%Backtracking linesearch parameters
max_bk_it = pars.max_bkt;
eta       = pars.eta_bkt;
%Step length parameters
tau       = pars.tau;
%Iterative solver tolerance
rtol      = 1e-16;
%Variables of the search direction
x_d = zeros(n,1);
t_d = 0;
y_d = zeros(m,1);
z_d = zeros(n,1);
k_d = 0;

rho = pars.rho;

if(rho == -1)
    rho = (n+1)-sqrt(n+1); %Feasibility weight in the potential function
end

sigma = (n+1)^2/rho;
%Error logs
err_function_count = size(pars.errFcn,2);
err_logs = zeros(max_iter,err_function_count);
if(pars.save_iters)
	iters_log = {};
end

%Calculate the initial centrality measure
mu   = x0'*z0+t0*k0;
mu   = mu/(n+1);

%Calculate the initial residuals
p_r  = t0*b - A*x0;
d_r  = t0*c - A'*y0 - z0;
g_r  = c'*x0 - b'*y0 + k0;
comp = x0'*z0+k0*t0;

n_p_r = norm(p_r);
n_p_r_t = n_p_r/t0;
n_d_r = norm(d_r);
n_d_r_t = n_d_r/t0;
n_g_r = abs(g_r);
gap_t = n_g_r/t0;

merit   = rho*log(comp) - sum(log(x0))-sum(log(z0))-log(t0)-log(k0);
g_merit = rho/comp*[zeros(m,1);z0;k0;x0;t0] - [zeros(m,1);1./x0;1./t0;1./z0;1./k0];

%rPrint the header
fprintf('%3s,%5s,%5s,%5s,%5s,%5s,%5s,%6s,%6s,%6s,%6s| %s \n','Itn','im','rtol','ib','is','pinf','dinf','cinf','gap','step','tau','errs');

%First iterates
x = x0;
t = t0;
y = y0;
z = z0;
k = k0;

%Centering parameter
compound_itn = 0;

solver  = pars.direction; 
for itn = 1:max_iter

    %Form the 3 by 3 system    
    x_sol     = K\RHS;
    
    %---------------------------------------------------------
    %Calculate step size to the boundary edge
    %--------------------------------------------------------- 
    alph_xz = find_step_to_neigh_bdry() 
    x       = x+alph_xz*d_x;
    t       = t+alph_xz*d_t;
    y       = y+alph_xz*d_y;
    z       = z+alph_xz*d_z;
    k       = k+alph_xz*d_k;
 
    y = alph_xz*x_sol(1:m);
    x = alph_xz*x_sol(m+1:m+n);
    t = alph_xz*x_sol(m+n+1);
    z = alph_xz*x_sol(n+m+2:2*n+m+1);
    k = alph_xz*x_sol(2*n+m+2);

    %-------------------------------------------------------
    %Newton Centering step
    %------------------------------------------------------- 
    for n_iter = 1:max_n_iter
        RHS  = build_centering_rhs();
        x_sol= K\RHS;                 
        
        centrality = eval_centrality(x,t,y,z,k);
        if centrality < pars.centrality_tol
            break;
        end
    end
    %--------------------------------------------------------
    %End newton centering step
    %---------------------------------------------------------
    
    %---------------------------------------------------------
    %Deal with failures
    %--------------------------------------------------------- 
    if fail
        fprintf('Linesearch failed \n');
        plot_present_point_and_direction(xs,ts,ys,zs,ks,x_d,t_d,y_d,z_d,k_d,...
        10000,A,b,c,max_alph_xz,rho,local_descent,eta,merit_a);
        %keyboard 
    end
    
    %--------------------------------------------------------
    %Update the iterates
    %-------------------------------------------------------- 

    if~fail
        %Calculate the new centrality measure
        comp   = x'*z+t*k;
        mu     = comp/(n+1);
        p_r    = t*b - A*x;
        d_r    = t*c - A'*y - z;
        g_r    = c'*x-b'*y+k;
        n_p_r_t = norm(p_r)/t;
        n_d_r_t = norm(d_r)/t;
        gap_t   = (c'*x-b'*y+k)/t; 
		
        %save the iterate
		if(pars.save_iters)
				iters_log = {iters_log{:},[x;t;y;z;k]};
		end

        g_merit = rho/comp*[zeros(m,1);z;k;x;t] - [zeros(m,1);1./x;1./t;1./z;1./k];
    end
    

    %--------------------------------------------------------------
    % Print and log
    %---------------------------------------------------------------
    
    %---------Append err functions that require inner info 
    pars.errFcn{err_function_count+1} = @(s)(max(x.*z)/min(x.*z));
    pars.errFcn{err_function_count+2} = @(s)(local_descent);
    pars.errFcn{err_function_count+3} = @(s)(t);    
    %pars.errFcn{err_function_count+2} = @(x)(inner_product_bound);  
    %---------Log the errors and print---------- 
    s_log = '%3i,%5i,%6.1f,%5i,%3i,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f |';
        for j=1:size(pars.errFcn,2)
           %Evaluate the errors
           err_logs(itn,j) = pars.errFcn{j}([x;t;y;z;k]);
           s_log = [s_log,', %3.3d'];
         end

    s_log = [s_log,'\n'];
    if size(pars.errFcn,2) > 0
        fprintf(s_log,itn,itn_aff,log10(rtol),bk_it,...
                istop,log10(n_p_r_t),log10(n_d_r_t),...
                log10(comp),log10(gap_t),log10(alfa),log10(full(t_d)),err_logs(itn,:));
    else
        fprintf(s_log,itn,itn_aff,log10(rtol),bk_it,...
                istop,log10(n_p_r_t),log10(n_d_r_t),...
                log10(comp),log10(gap_t),log10(alfa),log10(full(t_d))); 
    end 
    %---------End of logging----------------------
    
    %---------Evaluate exit criteria-------------------
    if(compound_itn >= pars.max_compound)
      fprintf('Maximum compound iterations reached\n');
      break;
    end
    
    if fail
       break;
    end
    
    if(n_p_r_t < tol & n_d_r_t < tol & (n+1)*mu < tol & gap_t < tol)
        break;
    end
    if(~isempty(pars.stopFcn))
        if(pars.stopFcn([x;y;z],pars.errFcn))
            break;
        end
    end
    %-------end of stop criteria --------------------
end %------ End of main iteration

%--------------------------------------------
%Functions that need the workspace
%--------------------------------------------
function K = Build_System()
    %Define the primal hessian at x
    Hess_x = primal_hessian(x);
    %Define the handle for the 3 by 3 system
     K = [[sparse(m,m),A          ,-b         ,sparse(m,n)       ,sparse(m,1)];
          [-A'        ,sparse(n,n),c          ,-speye(n,n)       ,sparse(n,1)];...
          [b'         ,-c'        ,sparse(1,1),sparse(1,n)       ,-1          ];...
          [sparse(n,m),mu*Hess_x  ,sparse(n,1),speye(n,n)        ,sparse(n,1)];...
          [sparse(1,m),sparse(1,n),mu/t^2     ,sparse(1,n)       ,1         ]];
  
end

function RHS = build_predictor_rhs()
   %Define the RHS for the affine direction
    RHS = [p_r;d_r;g_r;-[z;k]];  
end

function RHS = build_centering_rhs()
   %Define the RHS for the affine direction
    RHS = [sparse(n,1);sparse(m,1);sparse(1);-[z;k]+[mu*Hess_x*x;mu/t^2]];  
end


%--------------------------------------------
%End of functions that need the prip workspace 

end %end of nscs
end
