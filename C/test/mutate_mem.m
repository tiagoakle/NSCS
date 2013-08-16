%Exercises the mutate mem function to test 

%Load the library
[ret,warn] = loadlibrary('../Lib/liblinear_solvers','../Include/linear_solvers.h');
%Add the path to the matlab function that makes the call
addpath '../src'
%Generate a random system 
n = 5;
A = [[2,3,0,0,0];[3,0,4,0,6];[,0,-1,-3,2,0];[0,0,1,0,0];[0,4,2,0,1]];
b = [8, 45, -3, 3, 19]';
b2 = [1,1,1,1,1]';

%We assume that the library with name liblinear_solvers has been loaded
[I,J,V] = find(A); %Convert to triplet form
I = I-1;           %Shift to c notation
J = J-1;
n       = size(A,1);
nnz     = size(V,1); %We are using a dense matrix

%Desperate test make everything a pointer
lpV     = libpointer('doublePtr',V);
lpI     = libpointer('int32Ptr',I);
lpJ     = libpointer('int32Ptr',J);
lpX     = libpointer('doublePtr',zeros(n,1));
lpB     = libpointer('doublePtr',b);
ret     = calllib('liblinear_solvers','mutate_mem',lpX,lpI,lpV,lpB,nnz,n);
x       = lpX.Value;
