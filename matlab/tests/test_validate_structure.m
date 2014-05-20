
%To write the errors to 
fid = fopen('./tmp.out','w');
msg = 'Error: test validate structure ';

%Test the validate problem function 
problem = struct;
problem.n_pos = 0;
problem.n_exp_cones = 0;
problem.n_free = 10;
problem.n = 10;
problem.m = 11;
problem.A = zeros(11,10);
problem.c = zeros(10,1);
problem.b = zeros(11,1);

msg2 = 'Check correct problem passes test';
[r,problem]=validate_problem_structure(problem,fid);
assert(r==0,[msg,msg2]);

msg2= 'Check A present';
problemA=rmfield(problem,'A');
[r,problem]=validate_problem_structure(problemA,fid);
assert(r==-8,[msg,msg2]);

msg2= 'Check A size';
problem.A = zeros(10,10);
[r,problem]=validate_problem_structure(problem,fid);
assert(r==-8,[msg,msg2]);
problem.A = zeros(11,10);

msg2= 'Check b present';
problemb = rmfield(problem,'b');
[r,problem]=validate_problem_structure(problemb,fid);
assert(r==-8,[msg,msg2]);

msg2= 'Check b size';
problem.b = zeros(10,1);
[r,problem]=validate_problem_structure(problem,fid);
problem.b = zeros(11,1);
assert(r==-8,[msg,msg2]);

msg2= 'Check c present';
problemc=rmfield(problem,'c');
[r,problem]=validate_problem_structure(problemc,fid);
assert(r==-8,[msg,msg2]);

msg2= 'Check c size';
problem.c = zeros(11,1);
[r,problem]=validate_problem_structure(problem,fid);
assert(r==-8,[msg,msg2]);
problem.c = zeros(10,1);


fclose(fid);
