%script to test the build kkt routine
%load the library
[ret,warn] = loadlibrary('../lib/liblinear_solvers','../include/linear_solvers.h');
%add the path to the matlab function that makes the call
addpath '../src'

load './test_data/minentropy_system.mat'
[m,n] = size(A);
%%Select the size of the problem to generate
%m = 10;
%n = 100;
%density = 0.7;
%%Generate a constraint matrix
%A = sprand(m,n,density);
%L = sparse(n,n);
%%Generate a lower bidiagonal matrix L and form H = LL'
%L = diag(sparse(rand(n,1)))+diag(sparse(rand(n-1,1)),-1);
%H = L*L';

%Choose the regularization paramters and centrality value
mu = 0.5;
gamma = 1.e-10;
delta = 1.e-10;

[aI,aJ,aV] = find(A);
nnzA       = size(aI,1);
[hI,hJ,hV] = find(H);
nnzH       = size(hI,1);
finalnnz   = m + 2*nnzA+nnzH;

lpI     = libpointer('int32Ptr',zeros(1,finalnnz));
lpJ     = libpointer('int32Ptr',zeros(1,finalnnz));
lpV     = libpointer('doublePtr',zeros(1,finalnnz));

%int form_kkt_system_triplets(double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, int* aV, int nnzA, int m, int n, double delta, double gamma, int* I, int* J, double* V);

ret     = calllib('liblinear_solvers','form_kkt_system_triplets',mu,hI-1,hJ-1,hV,nnzH,aI-1,aJ-1,aV,nnzA,m,n,delta,gamma,lpI,lpJ,lpV);
KKT     = sparse(double(lpI.Value)+1,double(lpJ.Value)+1,lpV.Value);
KKT2    = [[delta*speye(m),A];[A',-mu*H-gamma*speye(n)]]; 

fprintf('Max difference %d\n', full(max(max(abs(KKT-KKT2)))));

%NOW test the function that returns the triplets in cco
lppI     = libpointer('int32PtrPtr',1);
lppP     = libpointer('int32PtrPtr',1);
lppV     = libpointer('doublePtrPtr',1.0);
ret     = calllib('liblinear_solvers','form_kkt_system',mu,hI-1,hJ-1,hV,nnzH,aI-1,aJ-1,aV,nnzA,m,n,delta,gamma,lppI,lppP,lppV);
fprintf('form_kkt_system returned %i\n',ret);

%------- Test 3 --------------------------------------
%Factor the matrix
%int factor_kkt_system(void** Numeric, int** Ai, int** Ap, double** Av, int n)
numeric  = libpointer('voidPtrPtr');
ret     = calllib('liblinear_solvers','factor_kkt_system',numeric,lppI,lppP,lppV,n+m);
fprintf('factor_kkt_system returned %i\n',ret);
    
%------- Part 4 ----------------------------------------
%Call the routine to solve 
%int solve_factored_system(void** Numeric, int** Ai, int** Ap, double** Av, int m, int n, int nnz, double* rhs, double* x);
lpX       = libpointer('doublePtr',zeros(n+m,1));
rhs       = zeros(n+m,1); %Make a random rhs
rhs(1:m)  = b;
rhs(m+1:end) = -c;
ret       = calllib('liblinear_solvers','solve_factored_system',numeric,lppI,lppP,lppV,rhs,lpX);
fprintf('solve_factored_system returned %i\n',ret);
%compare the solution to the matlab solution

x_2    = KKT2\rhs;
fprintf('Solution difference %g\n',norm(lpX.value-x_2));

%------- Part 5 --------------------------------------
%Test the call to free the memory

%int free_factorization(void** Numeric, int** Ai, int** Ap, double** Av)
ret     = calllib('liblinear_solvers','free_factorization',numeric,lppI,lppP,lppV);

unloadlibrary liblinear_solvers
    
