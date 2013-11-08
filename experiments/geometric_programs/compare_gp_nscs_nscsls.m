%Runs a set of GPs in coneopt, nscs and nscs long_step
clear all
%Runs all GP problems in coneopt
%Saves the results to the table stored in 
% the tex file 'coneopt_gp_run_results.tex'
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
'jha88.eo'};
%'mra01.eo'}; %Commented out for speed

%Cell array for the results
results = {{'Prob name','KKT coneopt','Coneopt Status','nscs ls KKT','nscs ls status','nscs kkt','nscs status'}};

problem_count = size(problem_names,2);
for(j =1:problem_count)
    
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
    pars.n = 3*num_ter + 2*num_var + 1;
    pars.m = 2*num_ter + num_con+1;
    pars.echo = 4;

%    pars.secord = 1;
    pars.cnbfgsstps = 3;
    pars.theta = 0.7;
    pars.eta   = 0.5;
    pars.beta  = 0.2;
   
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
   
    %Algorithm parameters
    pars.solve_second_order = true;

    x0c  = [t00;up0;um0;w00;v00;y00];
    x0f  = [];

    nscs 
    nscs_kkt = state.kkt_solves; 
    nscs_sta = state.exit_reason;
     
    %--------------------------------------------------------------------------
    % Solve with nscs long step
    %--------------------------------------------------------------------------
    nscs_long_step
    nscs_ls_kkt = state.kkt_solves; 
    nscs_ls_sta = state.exit_reason;
     
    %--------------------------------------------------------------------------
    % Save the results 
    %--------------------------------------------------------------------------
    problem_result = {problem_names{j},cone_kkt,cone_sta,...
                      nscs_ls_kkt,nscs_ls_sta,...
                      nscs_kkt,nscs_sta};
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
clear -REGEXP '^(?!.*?results).*'
%Save to file
save 'compare_gp_run_results'

