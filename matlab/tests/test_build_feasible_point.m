%Test build feasible point 
msg = 'Err: build feasible point ';

  initial_ex=[-0.827838387734678;...
                     1.290927717327200;...
                     0.805102001257750];       
problem = struct;
problem.n_free = 0;
problem.n_pos  = 0;
problem.n_exp_cones = 0;

%Build an empty problem 
x0c = build_feasible_primal_point(problem);
assert(isempty(x0c),msg);

%Build a problem with free vars only
problem.n_free = 1;
x0c = build_feasible_primal_point(problem);
assert(x0c(1)==0,msg);
assert(length(x0c)==1,msg);

%Build one with free and non negative 
problem.n_free = 1;
problem.n_pos  = 1;
x0c = build_feasible_primal_point(problem);
assert(x0c(1)==0,msg);
assert(x0c(2)==1,msg);
assert(length(x0c)==2,msg);

%Build one with non negative only
problem.n_free = 0;
problem.n_pos  = 1;
x0c = build_feasible_primal_point(problem);
assert(x0c(1)==1,msg);
assert(length(x0c)==1,msg);

%Build one with only exponential cones
problem.n_free = 0;
problem.n_pos  = 0;
problem.n_exp_cones = 1;
x0c = build_feasible_primal_point(problem);
assert(all(x0c==initial_ex),msg);
assert(length(x0c)==3,msg);

%Build one with pos and exponential cones
problem.n_free = 0;
problem.n_pos  = 1;
problem.n_exp_cones = 1;
x0c = build_feasible_primal_point(problem);
assert(x0c(1)==1,msg);
assert(all(x0c(2:end)==initial_ex),msg);
assert(length(x0c)==4,msg);

%Build one with all 3 types
problem.n_free = 1;
problem.n_pos  = 1;
problem.n_exp_cones = 1;
x0c = build_feasible_primal_point(problem);
assert(x0c(1)==0,msg);
assert(x0c(2)==1,msg);
assert(all(x0c(3:end)==initial_ex),msg);
assert(length(x0c)==5,msg);

msg2 = 'Make sure the order of the exp cones is correct';
problem.n_free = 0;
problem.n_pos  = 0;
problem.n_exp_cones = 2;
x0c = build_feasible_primal_point(problem);
assert(all(x0c(1:2:end)==initial_ex),[msg,msg2]);
assert(all(x0c(2:2:end)==initial_ex),[msg,msg2]);
assert(length(x0c)==6,msg);


