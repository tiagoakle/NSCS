%Script to test linear solvers
%Load the library
[ret,warn] = loadlibrary('../Lib/liblinear_solvers','../Include/linear_solvers.h');
%Add the path to the matlab function that makes the call
addpath '../interface'
%Generate a random system 
n = 5;
A = [[2,3,0,0,0];[3,0,4,0,6];[,0,-1,-3,2,0];[0,0,1,0,0];[0,4,2,0,1]];
b = [8, 45, -3, 3, 19]';
b2 = [1,1,1,1,1]';
%Call the solver 
x = linear_solver(A,b);
x2 = linear_solver(A,b2);
%Release the library
unloadlibrary liblinear_solvers
