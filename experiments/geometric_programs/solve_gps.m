
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


%NSCS parameters

%Generates the parameters structure for nscs long step with the default parameters
nscs_pars = set_default_pars_nscs_long_step();
%Set up some parameters
   nscs_pars.stop_primal= 1e-7;                 %Stopping criteria p_res/rel_p_res<stop_primal.
   nscs_pars.stop_dual  = 1e-7;
   nscs_pars.stop_gap   = 1e-6;
   nscs_pars.stop_mu    = 1e-7;
   nscs_pars.stop_tau_kappa = 1.e-5;
  
%Cell array for the results
results = {{'Prob name','KKT coneopt','Coneopt Status','nscs lsnt kkt','nscs lsnt status','msk iter','msk status'}};
problem_count = size(problem_names,2);
for(j =1:problem_count)
    fprintf('Will solve problem %s \n',problem_names{j});
    problem_file_name = problem_names{j};
    %Add the path to the file
    problem_file_name = ['./gp/',problem_file_name];
    [AA,bb,cc,num_ter,num_var,num_con,mskA,mskConstraints] = read_gp(problem_file_name);
    
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
   
    %x0c  = [t00;up0;um0;w00;v00;y00];

    x0c  = [];
    x0f  = [];

    %--------------------------------------------------------------------------
    % Solve with nscs long step with nt scaling
    %--------------------------------------------------------------------------
    pars.use_nesterov_todd_scaling = true;
    [xc,xf,y,z,t,k,info] = nscs_long_step(problem,x0f,x0c,nscs_pars);
    nscs_lsnt_kkt = info.kkt_solves; 
    nscs_lsnt_sta = info.exit_reason;

    %Evaluate the solution on the original probem
    nscs_x = xc/t;
    nscs_x = nscs_x(2:num_var+1)-nscs_x(num_var+2:2*num_var+1);
    nscs_ter = mskA*nscs_x;
    nscs_eu  = exp(nscs_ter-bb(1:num_ter));
    const_vals = mskConstraints*nscs_eu;
    nscs_obj    = const_vals(1)
    nscs_const  = const_vals(2:end)

    %--------------------------------------------------------------------------
    % Save the results 
    %--------------------------------------------------------------------------
    problem_result = {problem_names{j},...
                      nscs_lsnt_kkt,nscs_lsnt_sta};
    results = {results{:},problem_result};
 

end


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



