clear all
%Problem 'fiac81b.eo' is detected as infeasible by nscs while optimal by mosek.
%Evaluate if the resulting direction is in fact a recession direction!

%Check for linear feasibility and dual conic feasibility.

problem_file_name = 'fiac81b.eo'; 

%Build the conic version of the problem
problem_file_path = ['./gp/',problem_file_name];
[AA,bb,cc,num_ter,num_var,num_con] = read_gp(problem_file_path);
 

%--------------------------------------------------------------------------------
%NSCS parameters

%Generates the parameters structure for nscs long step with the default parameters
nscs_pars = struct;
%Set up the default parameters
   nscs_pars.max_iter   = 1000;  %Maximum outer iterations
   nscs_pars.max_affine_backtrack_iter = 300;    %Maximum affine backtracking steps
   nscs_pars.backtrack_affine_constant = 0.9;   %Affine backtracking constant

    %XXX: changed from 0.98 for gp testing
   nscs_pars.eta        = 0.8;                %Multiple of step to the boundary
   nscs_pars.stop_primal= 1e-7;                 %Stopping criteria p_res/rel_p_res<stop_primal.
   nscs_pars.stop_dual  = 1e-7;
   nscs_pars.stop_gap   = 1e-7;
   nscs_pars.stop_mu    = 1e-12;
   nscs_pars.stop_tau_kappa = 1.e-5;
   nscs_pars.solve_second_order = false;

   nscs_pars.print      = 2;                     %Level of verbosity from 0 to 11
   %Regularization for the linear solver
   nscs_pars.delta      = 5e-5;
   nscs_pars.gamma      = 5e-5;
   nscs_pars.max_iter_ref_rounds = 20;

%--------------------------------------------------------------------------------
% Set mosek parameters

%Mosek paths and parameters
mosek_path = '~/mosek/7/tools/platform/osx64x86/bin/mskexpopt';
mosek_parameter_file = '/mosek_gp_pars'
current_dir = cd;
problem_base_path = [current_dir,'/gp/'];
parameter_file_path = [current_dir,mosek_parameter_file];


    %--------------------------------------------------------------------------
    % Solve with nscs long step with nt scaling
    %--------------------------------------------------------------------------
    %Set up call to nscs
    % starting point:
    t00 = 1;
    up0 = ones(num_var,1);
    um0 = ones(num_var,1);
    w00 = -ones(num_ter,1);
    v00 = ones(num_ter,1);  
    y00 = 0.5*ones(num_ter,1); 
    x0c  = [t00;up0;um0;w00;v00;y00];
    x0f  = [];


    pars.use_nesterov_todd_scaling = true;
    problem = struct;
    problem.A = AA;
    problem.b = bb;
    problem.c = cc;
    problem.n = 2*num_var + 3*num_ter +1;
    problem.m =2*num_ter + num_con+1;
    problem.n_free = 0;
    problem.n_constrained = 2*num_var+1+3*num_ter;
    problem.n_pos       = 2*num_var+1;
    problem.n_exp_cones   = num_ter;
    

    [xc,xf,y,z,t,k,info] = nscs_long_step(problem,x0f,x0c,nscs_pars);
    nscs_lsnt_kkt = info.kkt_solves; 
    nscs_lsnt_sta = info.exit_reason;
         
    %-------------------------------------------------------------------------
    % Solve with MOSEK
    %-------------------------------------------------------------------------
    %Call mosek on each of the problems
 
    problem_path = [problem_base_path,problem_file_name];
    system_call  = [mosek_path,' -p ',parameter_file_path,' ',problem_path];
    %System call to mosek solver
    [r,s] = system(system_call);
    %S contains the sysout output
    %Extract the iteration count (find all lines that start with a number and extract the last) 
    it_lines = regexp(s,'\n\d');
    it_count_ix = it_lines(end);
    mosek_iteration_count = str2num(s(it_count_ix:it_count_ix+2))
    %Extract the status
    [tokens,t] = regexp(s,'Interior-point solution summary\n  Problem status  : (?<problemstatus>\S*)\n  Solution status : (?<solstat>\S*)','tokens');
    %Extract the status of the dual solution
    msk_feas = tokens{1}{1};
    mosek_status = tokens{1}{2};
   

%Check the dual feasibility of z.  
% z is formed by first 2*num_var variables for the free variables, then 3*num_ter for the exponential cones
% We build a problem struct and use it 
eval_dual_feas(problem,z)
