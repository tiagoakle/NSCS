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
  for problem_index = 4
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

    [M,N]   = size(P.A);
    b       = P.b(:);
    M2      = length(P.b);
    
    if M ~= M2
        error('A and b do not match');
    end
    
    pars.n = 3*N;
    pars.m = N+M;
    
    % build cone:
    K.npos = 0;
    K.npow = 0;
    K.nexp = N;
    K.nlog = 0;
    K      = getbarrpar(K);
    
    % build A:
    nnzA    = nnz(P.A);
    AA      = sparse([],[],[],pars.m,pars.n,nnzA+N);
    tmp     = N*(N+M);
    tmp     = tmp + 1:(N+M+1):(2*N*(N+M)-M);
    AA(tmp) = 1;
    tmp     = N+1:pars.m;
    tmp2    = 2*N+1:pars.n;
    AA(tmp,tmp2) = P.A;
    
    % build b:
    bb             = ones(pars.m,1);
    bb(N+1:pars.m) = b;
    bb             = sparse(bb);
    
    % build c:
    cc      = zeros(pars.n,1);
    cc(1:N) = -ones(N,1);
    cc      = sparse(cc);
 
  %Problem parameters
  
  fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,m,n);
 
 %  %-----------------------------
 %  % Call coneopt 
 %  %-----------------------------
 %  K = struct;
 %  % build cone:
 %  K.npos = 0;
 %  K.npow = 0;
 %  K.nexp = n;
 %  K.nlog = 0;
 %  K      = getbarrpar(K);
 %  
 %  %parameters 
 %  pars   = struct;
 %  pars.n = size(AA,2);
 %  pars.m = size(AA,1);
   pars.echo = 4;
   pars.secord = 0;
   pars.centmeastype = 5;
    
   pfeas = [-1.051383945322714;
              1.258967884768947;
              0.556409619469370];

     
    % starting point:
    u0  = -ones(n,1);
    v00 = ones(n,1);  
    x0  = 0.5*ones(n,1);
    %initial point 
    v0 = struct;
    v0.x  = [u0;v00;x0];
    
    pars.rhoP = 1.e-6;
    pars.rhoD = 1.e-6;
    pars.rhoA = 1.e-7;
  
    pars.echo   = 4;
    pars.beta   = 0.99;
    pars.trace  = 3;
    pars.secord = 1; 
   
    pars.permuteM = 0;
    % call to coneopt:
    R = coneopt(AA,bb,c,v0,K,pars);
 
   results = {results{:},{P.name,R.dat.nkktsolves,R.status}}; 
   %Prepare the call to nscs long step 
   %Extract the problem data and build the problem structure
   problem = struct;
   problem.A = AA;
   problem.b = bb;
   problem.c = cc;
   
   %Problem parameters
   problem.m = size(AA,1);
   problem.n = size(AA,2);
   problem.n_free = 0;
   problem.n_pos = 0;
   problem.soc_cones   = 0;
   problem.n_soc_cones = 0;
   problem.n_sdp_cones = 0;
   problem.n_exp_cones   = n;
   
   pars = set_default_pars_nscs_long_step(); 
   pars.solve_second_order = false;
   pars.neigh = 0.99;

   pars.stop_primal   = 1e-6;                 
   pars.stop_dual     = 1e-6;
   pars.stop_gap_res  = 1e-6; 
   pars.stop_gap      = 1e-7;
   pars.stop_tau_kappa= 1e-7;

    pars.stop_primal   = 1e-7;                 %stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual     = 1e-7;
    pars.stop_gap_res  = 1e-7;
    pars.stop_gap   = 1e-8;
    pars.stop_tau_kappa = 1.e-7;

   pars.print = 1;
   x0f        = []; 
   [xc,xf,y,s,t,k,info] = nscs_long_step(problem,x0f,v0.x,pars);
 end 

res = results{1};
fprintf('%s %s %s\n',res{1},res{2},res{3});
for(j=2:length(results))
    res = results{j};
    fprintf('%s %i %s\n',res{1},res{2},res{3});
end



