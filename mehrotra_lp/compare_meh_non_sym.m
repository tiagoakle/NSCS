%Loads an lp from lpnetlib and then 
% calls the mehrotra solver and the non symmetric long step
clear all
  
  addpath '../matlab'
  addpath '../coneopt'
  %Load the file that contains the indices for the
  %ufget netlib lps which are in standard form
  load 'standard_form_indices.mat' 

  %Choose a problem from the list
  problem_index = 7;
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
  A = P.A;
  b = P.b;
  c = P.aux.c;
fprintf('Problem %s \n',prob_name); 
[x,y,s,info] = mehrotra_lp_solver(A,b,c);
[x_ns,y_ns,s_ns,info_none] = non_symmetric_long_step(A,b,c,'linear');
[x_ns,y_ns,s_ns,info_ns] = non_symmetric_long_step(A,b,c,'arc_search');
[x_ns,y_ns,s_ns,info_ns] = non_symmetric_long_step(A,b,c,'secord');
%
%%now call ecos 
%n = size(A,2);
%G = -speye(n);
%h = zeros(n,1);
%dims.l = n;
%dims.q = [];
%[x_ecos,y_ecos,info_ecos,s_ecos,z_ecos] = ecos(-c,G,h,dims,A,b);
%
%%Now call sedumi
%[x_s,y_s] = sedumi(A,b,c);
