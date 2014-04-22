%Runs a problem in standard form from lp netlib 
% in the potential reduction hsd primal dual.. 

clear all

%Load the file that contains the indices for the
%ufget netlib lps which are in standard form
load 'standard_form_indices.mat' 

problem_count = length(st_ix);
problem_ix = 2;

%Get the ufget id 
problem_uf_ix = st_ix(problem_ix);
%Get the problem from ufget
P = UFget(problem_uf_ix);  
%Extract the name
prob_name = [P.name];
%Substitute front slash for space
prob_name(find(prob_name=='/'))=' ';

%Extract the problem data and build the problem structure
A = P.A;
b = P.b;
c = P.aux.c;
[m,n] = size(A);

%Call the solver
[x,y,s,info] = potential_reduction_neg_s(A,b,c); 
