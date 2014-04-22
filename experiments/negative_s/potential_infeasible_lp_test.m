%This test generates an unbounded LP and runs
%the algorithm that allows the dual slacks to become negative. 

clear all
close all

n = 100; %Size of the cone
m = 50; %Number of costraints

%-------
%Generate a feasible recession direction
x = rand(n,1);
x(n) = 1; %This will make sense after the following lines

%Generate a random A
A = sprandn(m,n,0.1);
%Calculate a modification of the last column to make Ax=0
d = -A(:,1:n-1)*x(1:n-1); 
A(:,n) = d;

%Find c s.t. c'*x<1.e-4
c = randn(n,1);
ctx = c'*x;
while(abs(ctx)<1.e-4)
    c = randn(n,1);
    ctx = c'*x;
end

if(ctx > 0)
    c = -c;
end
    
%Calculate a b in the range of A
b = A*randn(n,1);


fprintf('Generated a random feasible LP with %i constraints and %i variables\n',m,n);
fprintf('Ax %f \n',norm(A*x));

clear x s

[x,y,s,tau,kappa,info] = potential_reduction_neg_s(A,b,c);


