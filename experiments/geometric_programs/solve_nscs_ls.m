%Solves one problem with nscs_ls

%Runs a set of GPs in coneopt, nscs and nscs long_step
clear all

%Runs a short selection of fast to solve problems 
% uses coneopt, nscs_long_step with nt scaling and 
% nscs_long_step with no nt scaling

%Saves the results to the table stored in 
% the mat file 'compare_gp_run_result.mat'
addpath '../../coneopt'
addpath '../../matlab'
%Short and easy
problem_name = 'beck751.eo';

%Cell array for the results
results = {{'Prob name','nscs lsnt kkt','nscs lsnt status','nscs lsnt neg s kkt','nscs lsnt neg sstatus'}};
    
    fprintf('Will solve problem %s \n',problem_name);
    problem_file_name = problem_name;
    %Add the path to the file
    problem_file_name = ['./gp/',problem_file_name];
    [AA,bb,cc,num_ter,num_var,num_con] = read_gp(problem_file_name);
    
    %--------------------------------------------------------------------------
    % Solve with nscs
    %--------------------------------------------------------------------------
    %Set up call to nscs
    % starting point:
    t00 = 1;
    up0 = ones(num_var,1);
    um0 = ones(num_var,1);
    w00 = -ones(num_ter,1);
    v00 = ones(num_ter,1);  
    y00 = 0.5*ones(num_ter,1);
    
    %Extract the problem data and build the problem structure
    problem = struct;
    problem.A = AA;
    problem.b = bb;
    problem.c = cc;
    
    %Problem parameters
    problem.m =2*num_ter + num_con+1;
    problem.n = 3*num_ter + 2*num_var + 1;
    problem.n_free = 0;
    problem.n_constrained = 2*num_var+1+3*num_ter;
    problem.n_pos       = 2*num_var+1;
    problem.soc_cones   = 0;
    problem.n_soc_cones = 0;
    problem.n_sdp_cones = 0;
    problem.sdp_cones     = 0;
    problem.n_exp_cones   = num_ter;
    problem.n_power_cones = 0;
   
    x0c  = [t00;up0;um0;w00;v00;y00];
    x0f  = [];

    %Populate the default structure 
    set_default_pars_nscs_long_step;
    pars.eta = 0.98; 
    pars.backtrack_affine_constant = 0.9;   %Affine backtracking constant
    pars.stop_primal= 1e-7;                 %Stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual  = 1e-7;
    pars.stop_gap   = 1e-7;
    pars.stop_mu    = 1e-7;
    pars.stop_tau_kappa = 1.e-7;

    %Regularization for the linear solver
    pars.delta      = 5e-5;
    pars.gamma      = 5e-5;
    pars.max_iter_ref_rounds = 10;


    %--------------------------------------------------------------------------
    % Solve with nscs long step with nt scaling
    %--------------------------------------------------------------------------
    pars.use_nesterov_todd_scaling = true;
    [xc,xf,y,z,t,k,info] = nscs_long_step(problem,x0f,x0c,pars);
    nscs_lsnt_kkt = info.kkt_solves; 
    nscs_lsnt_sta = info.exit_reason;
    
  %  %--------------------------------------------------------------------------
  %  % Solve with nscs long step with nt scaling and negative s
  %  %--------------------------------------------------------------------------
  %  pars.use_nesterov_todd_scaling = true;
  %  pars.print = 2
  %  [xc_ns,xf_ns,y_ns,z_ns,t_ns,k_ns,info] = nscs_long_step_neg_s(problem,x0f,x0c,pars);
  %  nscs_lsnt_kkt_ns = info.kkt_solves; 
  %  nscs_lsnt_sta_ns = info.exit_reason;
  %   
  %  %--------------------------------------------------------------------------
  %  % Save the results 
  %  %--------------------------------------------------------------------------
  %  j = 2;
  %  problem_result = {problem_name,...
  %                    nscs_lsnt_kkt,nscs_lsnt_sta,nscs_lsnt_kkt_ns,nscs_lsnt_sta_ns};
  %  results = {results{:},problem_result};
 

  %  %Print out
  %  fid = 1;
  %  
  %     fprintf(fid,'%s, %s, %s \n',results{1}{1},results{1}{2},...
  %                                                results{1}{3},results{1}{4},...
  %                                                results{1}{5});
  %  for j=2:size(results,2)
  %     fprintf(fid,'%10s, %3i, %15s \n',...
  %                                                results{j}{1},results{j}{2},...
  %                                                results{j}{3},results{j}{4},...
  %                                                results{j}{5});
  %end

