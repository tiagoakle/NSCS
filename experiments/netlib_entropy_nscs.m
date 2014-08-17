
  %Solves the entropy minimization problem 
  % min sum x_i log x_i s.t. Ax = b
  
  clear all
  addpath '../'
  %Load the file that contains the indices for the
  %ufget netlib lps which are in standard form
  load 'standard_form_indices.mat' 

  %Choose a problem from the list
  problem_index = 4;
  %Extract the problem 
  problem_uf_ix = st_ix(problem_index);
  %Get the problem from ufget
  P = UFget(problem_uf_ix);
    
  %Extract the name
  prob_name = [P.name];
  %Substitute front slash for space
  prob_name(find(prob_name=='/'))=' ';
  
  %Problem parameters
  m = size(P.A,1);
  n = size(P.A,2);
 
  fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,m,n);
  pars = set_default_pars_nscs_long_step(); 
  pars.solve_second_order = true;
  
  pars.print  = 2;

  x0f         = [];
  [xc] = minentropy_nscs(P.A,P.b,ones(size(P.A,2),1),pars);


