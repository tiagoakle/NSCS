clear all
addpath '../../matlab'
%This problem should be divergent, lets see what happens
%minimize exp(a'x) 
%st exp(a'x)<1
%The objective shoud be zero but the minimizer non finite
a = 1
m = 3
n = 6

A = [[a -a 0 -1 0 0];[0 0 1 0 1 0];[0 0 0 0 0 1]];
c = [0 0 0 0 1 0]';
b = [0;1;1]
e = exp(1)
x0 = [1,1,0.2,1/(2*e),0.75,1/(2*e)]';
%What will nscs do?

%Generates the parameters structure for nscs long step with the default parameters
nscs_pars = struct;
%Set up the default parameters
   nscs_pars.max_iter   = 100;  %Maximum outer iterations
   nscs_pars.max_affine_backtrack_iter = 300;    %Maximum affine backtracking steps
   nscs_pars.backtrack_affine_constant = 0.9;   %Affine backtracking constant

    %XXX: changed from 0.98 for gp testing
   nscs_pars.eta        = 0.98;                %Multiple of step to the boundary
   nscs_pars.stop_primal= 1e-12;                 %Stopping criteria p_res/rel_p_res<stop_primal.
   nscs_pars.stop_dual  = 1e-12;
   nscs_pars.stop_gap   = 1e-12;
   nscs_pars.stop_mu    = 1e-12;
   nscs_pars.stop_tau_kappa = 1.e-5;
   nscs_pars.solve_second_order = true;

   nscs_pars.print      = 1;                     %Level of verbosity from 0 to 11
   %Regularization for the linear solver
   nscs_pars.delta      = 5e-5;
   nscs_pars.gamma      = 5e-5;
   nscs_pars.max_iter_ref_rounds = 20;
   


    %--------------------------------------------------------------------------
    % Solve with nscs
    %--------------------------------------------------------------------------
    %Set up call to nscs
    
    %Extract the problem data and build the problem structure
    problem = struct;
    problem.A = A;
    problem.b = b;
    problem.c = c;
    
    %Problem parameters
    problem.m = 3;
    problem.n = 6;
    problem.n_free = 0;
    problem.n_constrained = 6;
    problem.n_pos       = 3;
    problem.soc_cones   = 0;
    problem.n_soc_cones = 0;
    problem.n_sdp_cones = 0;
    problem.sdp_cones     = 0;
    problem.n_exp_cones   = 1;
    problem.n_power_cones = 0;
   
    x0f  = [];

    %--------------------------------------------------------------------------
    % Solve with nscs long step with nt scaling
    %--------------------------------------------------------------------------
    pars.use_nesterov_todd_scaling = true;
    [xc,xf,y,z,t,k,info] = nscs_long_step(problem,x0f,x0,nscs_pars);
    nscs_lsnt_kkt = info.kkt_solves; 
    nscs_lsnt_sta = info.exit_reason;

    
    %Build a mosek file 
    %One constraint two terms one variable 
    %c0,c1 = 1
    %First term is the objective 
    %second term the constraint
    msk_header = [1,1,2,1,1,0,1];
    msk_data   = [[0,0,a];[1,0,a]];
    %Build the string
    msk_str = '';
    for j =1:size(msk_header,2)
        msk_str = [msk_str,sprintf('%i \n',msk_header(j))];
    end
    for j =1:size(msk_data,1)
        msk_str = [msk_str,sprintf('%i %i %g\n',msk_data(j,:))];
    end
    f = fopen('./gp/tmp.eo','w')
    fprintf(f,'%s',msk_str)
    fclose(f) 

    %-------------------------------------------------------------------------
    % Solve with MOSEK
    %-------------------------------------------------------------------------
    %Call mosek on each of the problems
    %Mosek paths and parameters
    mosek_path = '~/mosek/7/tools/platform/osx64x86/bin/mskexpopt';
    mosek_parameter_file = '/mosek_gp_pars'
    current_dir = cd;
    problem_base_path = [current_dir,'/gp/'];
    parameter_file_path = [current_dir,mosek_parameter_file];

    problem_name = 'tmp.eo';
    problem_path = [problem_base_path,problem_name];
    system_call  = [mosek_path,' -p ',parameter_file_path,' ',problem_path]
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
   

 
