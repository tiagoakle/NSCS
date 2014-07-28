%Experiment to find the solution to 
% min sum xlogx
% subject to Ax=b

%Loads an lp from lpnetlib and then 
% coneopt

clear all
  addpath '../'
  addpath '../coneopt/'
  addpath '../matlab'
  %Load the file that contains the indices for the
  %ufget netlib lps which are in standard form
  load 'standard_form_indices.mat'; 
  problem_count = length(st_ix);

  results = {{'Problem name','Iterations','Exit Flag'}}


  %Choose a problem from the list
  for problem_index = 1:5
    %Extract the problem 
    problem_uf_ix = st_ix(problem_index);
    %Get the problem from ufget
    P = UFget(problem_uf_ix);
      
    %Extract the name
    prob_name = [P.name];
    %Substitute front slash for space
    prob_name(find(prob_name=='/'))=' ';
    
    %Extract the problem data and build the problem structure
    %Problem parameters
    [m,n] = size(P.A);
    fprintf('%i: Loaded problem %s, m:%i, n:%i\n',problem_index,P.name,m,n);
  

    % min sum(xlog(x)) st Ax=b
    % using u,v,w and wexp(u/w)<v => u < wlog(v/w) = u<-wlog(w/v) => wlog(w/v)<-u
    % minimize -u 0 0  
    % create the constraints Aw=b, v = 1 u,v,w \in K_exp
    % [0 I 0] [u]   = 1
    % [0 0 A] [v]     b
    %         [w]
    
    AA =[ [sparse(n,n) speye(n,n)  sparse(n,n)];
          [sparse(m,n) sparse(m,n)        P.A ]];

    c  = [-ones(n,1);zeros(n,1);zeros(n,1)];
    bb = [ones(n,1);P.b];
 
  %Problem parameters
  
  fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,m,n);
 
   %-----------------------------
   % Call coneopt 
   %-----------------------------
   K = struct;
   % build cone:
   K.npos = 0;
   K.npow = 0;
   K.nexp = n;
   K.nlog = 0;
   K      = getbarrpar(K);
   %parameters 
   pars   = struct;
   pars.n = size(AA,2);
   pars.m = size(AA,1);
   pars.echo = 4;
   pars.secord = 0;
   pars.centmeastype = 5;
    
   pfeas = [-1.051383945322714;
              1.258967884768947;
              0.556409619469370];

    %Make the initial points
    x0c = [pfeas(1)*ones(n,1);pfeas(2)*ones(n,1);pfeas(3)*ones(n,1)];
    
   %initial point 
   v0 = struct;
   v0.x = x0c;
   pars.rhoP = 1.e-6;
   pars.rhoD = 1.e-6;
   pars.rhoA = 1.e-7;
   % call to coneopt:
   R = coneopt(AA,bb,c,v0,K,pars);
   results = {results{:},{P.name,R.dat.nkktsolves,R.status}}; 
   %Prepare the call to nscs long step 
   %Extract the problem data and build the problem structure
   problem = struct;
   problem.A = AA;
   problem.b = bb;
   problem.c = c;
   %Problem parameters
   problem.m = size(AA,1);
   problem.n = size(AA,2);
   problem.n_free = 0;
   problem.n_pos = 0;
   problem.soc_cones   = 0;
   problem.n_soc_cones = 0;
   problem.n_sdp_cones = 0;
   problem.n_exp_cones   = n;
   problem.print = 1;
   
   pars = set_default_pars_nscs_long_step(); 
   pars.second_order = false;
   pars.neigh = Inf;
    x0f        = [];
   %[xc,xf,y,s,t,k,info] = nscs_long_step(problem,x0f,x0c,pars);
 end 

res = results{1};
fprintf('%s %s %s\n',res{1},res{2},res{3});
for(j=2:length(results))
    res = results{j};
    fprintf('%s %i %s\n',res{1},res{2},res{3});
end



