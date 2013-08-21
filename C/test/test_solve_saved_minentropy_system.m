%Load the saved matrix 
load 'minentropy_system'
clear all
load minentropy_system
%load the library
[ret,warn] = loadlibrary('../lib/liblinear_solvers','../include/linear_solvers.h');
%add the path to the matlab function that makes the call
addpath '../src'

[m,n]= size(A);
%Choose the regularization paramters
gamma = 1.e-10;
delta = 1.e-10;

%Get the triplet form of the matrix
[aI,aJ,aV] = find(A);
nnzA       = size(aI,1);
[hI,hJ,hV] = find(H);
nnzH       = size(hI,1);
finalnnz   = m + 2*nnzA+nnzH;

%Define the c variables

lppI     = libpointer('int32Ptr');
lppP     = libpointer('int32Ptr');
lppV     = libpointer('doublePtr');
lpX     = libpointer('doublePtr',zeros(n+m,1));

rhs     = zeros(n+m,1);
rhs(1:m)     =  b;
rhs(m+1:m+n) = -c;
%rhs     = [b;-c];

ppnumeric  = libpointer('voidPtr');

%Call the solver
ret     = calllib('liblinear_solvers','form_kkt_system',mu,hI-1,hJ-1,hV,nnzH,aI-1,aJ-1,aV,nnzA,m,n,delta,gamma,lppI,lppP,lppV);
ret     = calllib('liblinear_solvers','factor_kkt_system',ppnumeric,lppI,lppP,lppV,n+m);
ret     = calllib('liblinear_solvers','solve_factored_system',ppnumeric,lppI,lppP,lppV,rhs,lpX);

KKT2    = [[delta*speye(m),A];[A',-mu*H-gamma*speye(n)]]; 

%Build the matrix in MATLAB and solve with backslash for comparisson
x_2    = KKT2\rhs;
fprintf('Solution difference %g\n',norm(lpX.value-x_2));

ret     = calllib('liblinear_solvers','free_factorization',ppnumeric,lppI,lppP,lppV);

%Test one go call
lpX2       = libpointer('doublePtr',zeros(n+m,1));
%int build_and_solve_linear_system( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, double* sol, double * rhs)
ret     = calllib('liblinear_solvers','build_and_solve_linear_system',mu,hI-1,hJ-1,hV,nnzH,aI-1,aJ-1,aV,nnzA,m,n,delta,gamma,lpX2,rhs);
fprintf('Solution difference single call solve %g\n',norm(lpX2.value-x_2));

%Test build outside and do one call
[kI,kJ,kV] = find(KKT2);
nnzK       = size(kI,1);
lpX3     = libpointer('doublePtr',zeros(n+m,1));
ret     = calllib('liblinear_solvers','solve_linear_system',lpX3,kI-1,kJ-1,kV,nnzK,rhs,n+m);
fprintf('Solution difference build outside and call %g\n',norm(lpX3.value-x_2));


%Save the correct solution to a csv
write_vector_to_csv('solution.csv',x_2);
write_matrix_to_csv('H.csv',H);
write_matrix_to_csv('A.csv',A);

unloadlibrary liblinear_solvers
    
