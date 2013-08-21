int solve_linear_system(double* x, int* I, int* J, double* V, int nnz, double *b ,int n);
//This function forms the matrix to factorize and stores it in triplet form in Ai, Aj, Av
int form_kkt_system_triplets(double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int* I, int* J, double* V);
//This function forms the regularized KKT matrix and stores it in compressed col form. Uses form_kkt_system_triplets
int form_kkt_system         ( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int** Ai, int** Ap, double** Av);
//This function receives the matrix to factor in cco form and retunrs
int factor_kkt_system       (void** Numeric, int* Ai, int* Ap, double* Av, int n);
//Solves the linear system using the factorization
int solve_factored_system(void* Numeric, int* Ai, int* Ap, double* Av, double* rhs, double* x);
//Clear the memory used for the factored matrix
int free_factorization(void* Numeric, int* Ai, int* Ap, double* Av);
//One call to execute all the solve cycle
int build_and_solve_linear_system( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, double* sol, double * rhs);
//Debugging 
int loop(double *x, double* b, int n);
