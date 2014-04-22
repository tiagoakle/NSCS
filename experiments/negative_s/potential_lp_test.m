clear all
close all


n = 100; %Size of the cone
m = 50; %Number of costraints

%-------
%Generate a feasible primal point
x = rand(n,1);

%Generate a feasible dual point
s = rand(n,1);

%Generate a feasible and bounded primal dual problem
A = sprandn(m,n,0.1);
b = A*x;
s = rand(n,1);
c = A'*ones(m,1) + s;
fprintf('Generated a random feasible LP with %i constraints and %i variables\n',m,n);
clear x s

[x,y,s,info] = potential_reduction_neg_s(A,b,c);


