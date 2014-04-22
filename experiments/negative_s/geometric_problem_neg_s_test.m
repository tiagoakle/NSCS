%Runs a set of GPs in the experimental negative S code
clear all

addpath '../../matlab'
addpath '../geometric_programs'

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

PROBLEM_IX = 7;
problem_name = problem_names{PROBLEM_IX}

    fprintf('Will solve problem %s \n',problem_name);
    problem_file_name = problem_name;
    %Add the path to the file
    problem_file_name = ['../geometric_programs/gp/',problem_file_name];
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
    pars.max_iter = 1000;
%    %--------------------------------------------------------------------------
%    % Solve with the experimental negative s code
%    %--------------------------------------------------------------------------
     [xc,xf,y,s,info] = potential_reduction_neg_s_exponential_cone(problem,x0f,x0c,pars);

