%Script to test linear solvers
%Load the library
[ret,warn] = loadlibrary('liblinear_solvers','linear_solvers.h');
%Generate a random system 
n = 10;
A = randn(n);
b = randn(n,1);
%Call the solver 
x = linear_solver(A,b);
%Release the library
unloadlibrary liblinear_solvers
