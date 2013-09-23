#include "spmat.h"

#ifdef __cplusplus
extern "C" {
#endif

int solve_linear_system(double* x, int* I, int* J, double* V, int nnz, double *b ,int n);
void form_kkt_system_triplets( double mu,\
                              int* hI,\
                              int* hJ,\
                              double* hV,\
                              int nnzH,\
                              int* aI,\
                              int* aJ,\
                              double* aV,\
                              int nnzA,\
                              int m,\
                              int n,\
                              double delta,\
                              double gamma,\
                              int* I,\
                              int* J,\
                              double* V);
int form_kkt_system( double mu,\
                     int* hI,\
                     int* hJ,\
                     double* hV,\
                     int nnzH,\
                     int* aI,\
                     int* aJ,\
                     double* aV,\
                     int nnzA,\
                     int m,\
                     int n,\
                     double delta,\
                     double gamma,\
                     int** Ai,\
                     int** Ap,\
                     double** Av);
int factor_kkt_system(void** Numeric, int* Ai, int* Ap, double* Av, int n, int m);
int solve_factored_system(void* Numeric, int* Ai, int* Ap, double* Av, double* rhs, double* x);
int free_factorization(void* Numeric, int* Ai, int* Ap, double* Av);
int loop(double *x, double* b, int n);
int solve_kkt_system_no_structs(int m, int n,
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
                                double* dk );
int solve_kkt_system(double mu,\
                     spmat H,\
                     spmat A,\
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
                     double* dk );
void dummy_copy(double* Out, double* In, int n);

#ifdef __cplusplus
}
#endif
