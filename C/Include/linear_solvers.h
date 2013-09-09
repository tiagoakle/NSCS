//#include "nscs_sp.h"
int solve_linear_system(double* x, int* I, int* J, double* V, int nnz, double *b ,int n);
//This function forms the matrix to factorize and stores it in triplet form in Ai, Aj, Av
void form_kkt_system_triplets(double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int* I, int* J, double* V);
//This function forms the regularized KKT matrix and stores it in compressed col form. Uses form_kkt_system_triplets
int form_kkt_system         ( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int** Ai, int** Ap, double** Av);
//This function receives the matrix to factor in cco form and retunrs
int factor_kkt_system       (void** Numeric, int* Ai, int* Ap, double* Av, int n, int m);
//Solves the linear system using the factorization
int solve_factored_system(void* Numeric, int* Ai, int* Ap, double* Av, double* rhs, double* x);
//Clear the memory used for the factored matrix
int free_factorization(void* Numeric, int* Ai, int* Ap, double* Av);
//One call to execute all the solve cycle
//int build_and_solve_linear_system( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, double* sol, double * rhs);
//Debugging 
int loop(double *x, double* b, int n);
//Wrapper to build_and_solve_linear_system with no structs to call from matlab
int solve_kkt_system_no_structs(int m, int n,\
                                double mu,\
                                int* Hi,\
                                int* Hj,\
                                double* Hv,\
                                int Hnnz,\
                                int* Ai,\
                                int* Aj,\
                                double* Av,\
                                int Annz,\
                                double* b,\
                                double* c,\
                                double tau,\
                                double kappa,\
                                double delta,\
                                double gamma,\
                                double* r1,\
                                double* r2,\
                                double r3,\
                                double* r4,\
                                double r5,\
                                double* dy,\
                                double* dx,\
                                double* dt,\
                                double* ds,\
                                double* dk);

void dummy_copy(double* Out, double* In, int n);
//Builds and solves the 5x5 system
//int solve_kkt_system(double mu,\
//                     spmat H,\
//                     spmat A,\
//                     vec b,\
//                     vec c,\
//                     double tau,\
//                     double kappa,\
//                     double delta,\
//                     double gamma,\
//                     vec r1,\
//                     vec r2,\
//                     double r3,\
//                     vec r4,\
//                     double r5,\
//                     vec dy,\
//                     vec dx,\
//                     double* dt,\
//                     vec ds,\
//                     double* dk);
