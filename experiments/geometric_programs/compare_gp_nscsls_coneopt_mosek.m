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
problem_names = {...
'beck751.eo',...
'beck752.eo',...
'beck753.eo',...
'bss2.eo',...
'car.eo',...
'demb761.eo',...
'demb762.eo',...
'demb763.eo',...
'demb781.eo',...
'fang88.eo',...
'fiac81a.eo',...
'fiac81b.eo',...
'gptest.eo',...
'jha88.eo',...
'mra01.eo'}; %Commented out for speed

%Mosek paths and parameters
mosek_path = '~/mosek/7/tools/platform/osx64x86/bin/mskexpopt';
mosek_parameter_file = '/mosek_gp_pars'
current_dir = cd;
problem_base_path = [current_dir,'/gp/'];
parameter_file_path = [current_dir,mosek_parameter_file];

%NSCS parameters

%Generates the parameters structure for nscs long step with the default parameters
nscs_pars = struct;
%Set up the default parameters
   nscs_pars.max_iter   = 100;  %Maximum outer iterations
   nscs_pars.max_affine_backtrack_iter = 300;    %Maximum affine backtracking steps
   nscs_pars.backtrack_affine_constant = 0.95;   %Affine backtracking constant

    %XXX: changed from 0.98 for gp testing
   nscs_pars.eta        = 0.98;                %Multiple of step to the boundary
   nscs_pars.stop_primal= 1e-7;                 %Stopping criteria p_res/rel_p_res<stop_primal.
   nscs_pars.stop_dual  = 1e-7;
   nscs_pars.stop_gap   = 1e-7;
   nscs_pars.stop_mu    = 1e-12;
   nscs_pars.stop_tau_kappa = 1.e-5;
   nscs_pars.solve_second_order = true;

   nscs_pars.print      = 1;                     %Level of verbosity from 0 to 11
   %Regularization for the linear solver
   nscs_pars.delta      = 5e-5;
   nscs_pars.gamma      = 5e-5;
   nscs_pars.max_iter_ref_rounds = 20;



%Cell array for the results
results = {{'Prob name','KKT coneopt','Coneopt Status','nscs lsnt kkt','nscs lsnt status','msk iter','msk status'}};
problem_count = size(problem_names,2);
%for(j =1:problem_count)
j = 12    
    fprintf('Will solve problem %s \n',problem_names{j});
    problem_file_name = problem_names{j};
    %Add the path to the file
    problem_file_name = ['./gp/',problem_file_name];
    [AA,bb,cc,num_ter,num_var,num_con] = read_gp(problem_file_name);
    
    %----------------------------------------------------
    % coneopt call 
    %---------------------------------------------------
    % build cone:
    K.npos = 2*num_var+1;
    K.npow = 0;
    K.nexp = num_ter;
    K.nlog = 0;
    K      = getbarrpar(K);
     
    %Set the parameters
    pars = struct
    pars.n = 3*num_ter + 2*num_var + 1;
    pars.m = 2*num_ter + num_con+1;
    pars.echo = 4;

%    pars.secord = 1;
    pars.cnbfgsstps = 3;
    pars.theta = 0.7;
    pars.eta   = 0.8;
    pars.beta  = 0.8;
   
    pars.rhoP  = 1e-5;
    pars.rhoD  = 1e-5;
    pars.rhoA  = 1e-5;
    pars.rhoG  = 1e-5;
    pars.rhoI  = 1e-7;
    pars.rhoM  = 1e-7;

    %pars.centmeastype = 5;
    
    % starting point:
    t00 = 1;
    up0 = ones(num_var,1);
    um0 = ones(num_var,1);
    w00 = -ones(num_ter,1);
    v00 = ones(num_ter,1);  
    y00 = 0.5*ones(num_ter,1);
    
    v0.x  = [t00;up0;um0;w00;v00;y00];
    
    % call to coneopt:
    R = coneopt(AA,bb,cc,v0,K,pars);
    cone_kkt = R.dat.nkktsolves;
    cone_sta = R.status;
   
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

%    %--------------------------------------------------------------------------
%    % Solve with nscs long step no nt
%    %--------------------------------------------------------------------------
%    pars.use_nesterov_todd_scaling = false;
%    [xc,xf,y,z,info] = nscs_long_step(problem,x0f,x0c,pars);
%    nscs_ls_kkt = info.kkt_solves; 
%    nscs_ls_sta = info.exit_reason;

    %--------------------------------------------------------------------------
    % Solve with nscs long step with nt scaling
    %--------------------------------------------------------------------------
    pars.use_nesterov_todd_scaling = true;
    [xc,xf,y,z,t,k,info] = nscs_long_step(problem,x0f,x0c,nscs_pars);
    nscs_lsnt_kkt = info.kkt_solves; 
    nscs_lsnt_sta = info.exit_reason;
         
    %-------------------------------------------------------------------------
    % Solve with MOSEK
    %-------------------------------------------------------------------------
    %Call mosek on each of the problems
    problem_name = problem_names{j};
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
   

 
    %--------------------------------------------------------------------------
    % Save the results 
    %--------------------------------------------------------------------------
    problem_result = {problem_names{j},cone_kkt,cone_sta,...
                      nscs_lsnt_kkt,nscs_lsnt_sta,mosek_iteration_count,mosek_status};
    results = {results{:},problem_result};
 

%end


%Print out
fid = 1;

   fprintf(fid,'%s, %s, %s, %s, %s, %s, %s\n',results{1}{1},results{1}{2},...
                                              results{1}{3},results{1}{4},...
                                              results{1}{5},results{1}{6},...
                                              results{1}{7});
for j=2:size(results,2)
   fprintf(fid,'%10s, %3i, %15s, %3i, %15s, %3i, %15s \n',...
                                              results{j}{1},results{j}{2},...
                                              results{j}{3},results{j}{4},...
                                              results{j}{5},results{j}{6},...
                                              results{j}{7});
end

%Clear all but results 
%clear -REGEXP '^(?!.*?results).*'
%Save to file
%save 'compare_gp_run_results'


