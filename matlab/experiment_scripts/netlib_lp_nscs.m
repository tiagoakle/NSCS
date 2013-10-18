  clear all
  addpath '../'
  %Load the file that contains the indices for the
  %ufget netlib lps which are in standard form
  load 'standard_form_indices.mat' 

  %Choose a problem from the list
  problem_index = 5;
  %Extract the problem 
  problem_uf_ix = st_ix(problem_index);
  %Get the problem from ufget
  P = UFget(problem_uf_ix);
    
  %Extract the name
  prob_name = [P.name];
  %Substitute front slash for space
  prob_name(find(prob_name=='/'))=' ';
  
  %Extract the problem data and build the problem structure
  problem = struct;
  problem.A = P.A;
  problem.b = P.b;
  problem.c = P.aux.c;
  %Problem parameters
  problem.m = size(problem.A,1);
  problem.n = size(problem.A,2);
  problem.n_free = 0;
  problem.n_constrained = problem.n;
  problem.n_pos       = problem.n;
  problem.soc_cones   = 0;
  problem.n_soc_cones = 0;
  problem.n_sdp_cones = 0;
  problem.sdp_cones     = 0;
  problem.n_exp_cones   = 0;
  problem.n_power_cones = 0;
  x0c        = ones(problem.n,1);
  x0f        = []; 
  fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,problem.m,problem.n);
  nscs 

