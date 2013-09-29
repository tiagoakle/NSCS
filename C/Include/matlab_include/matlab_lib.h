#ifndef H_MATLAB
#define H_MATLAB

#define csi int
//int linesearch_atd(state_t state ,parameters_t  pars, problem_t prob);
#ifdef __cplusplus
extern "C" {
#endif

/**
 * Wrapper to call from matlab
 */
int linesearch_cent_no_structs( int m, int n, double*x, double*y, double*s, double tau, double kappa, double* dx,\
                                                                                                    double* dy,\
                                                                                                    double* ds,\
                                                                                                    double  dtau,\
                                                                                                    double  dkappa,\
                                                                                                    double  lsccent,\
                                                                                                    double  eta,\
                                                                                                    double  theta,\
                                                                                                    double  mu,\
                                                                                                    int max_backtrack,\
                                                                                                    int k_count,\
                                                                                                    int* nK,\
                                                                                                    int* tK,\
                                                                                                    double nu,\
                                                                                                    csi nnzH,\
                                                                                                    double* a,\
                                                                                                    int* nbacktrack,
                                                                                                    double* objvalb);

int linesearch_atd_no_structs( int m, int n, double* x, double* y, double* s, double tau, double kappa, double* dx,\
                                                                                                    double* dy,\
                                                                                                    double* ds,\
                                                                                                    double  dtau,\
                                                                                                    double  dkappa,\
                                                                                                    double  lscaff,\
                                                                                                    double  eta,\
                                                                                                    double  theta,\
                                                                                                    int max_backtrack,\
                                                                                                    int k_count,\
                                                                                                    int* nK,\
                                                                                                    int* tK,\
                                                                                                    double nu,\
                                                                                                    int nnzH,\
                                                                                                    double* a,\
                                                                                                    int* nbacktrack);

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

//Builds and solves the 5x5 system
void eval_hess_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* HI, int* HJ, double* HV);
void primal_feas_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* feas);
void dual_feas_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* feas);
void eval_grad_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, double* g);
void eval_cent_meas_no_structs(csi k_count,csi* nK, int* tK, double delta,\
                               double* x, double* s, csi n, csi nnzH, double mua,\
                               double* psi, double * hpsi, double* centmeas);


void nscs_no_structs(csi k_count, csi* nK, int* tK,\
                     csi *AI, csi* AJ, double* AV, csi nnzA, csi m, csi n,\
                     double *c, double *b,\
                     double* y, double* x, double* s, double* tau, double* kappa,\
                     int wy, int ws, int wt, int wk,\
                     int max_iter, int max_center_iter, double theta, double lscaff, double lsccent,\
                     double eta, int max_backtrack, double delta, double gamma, double beta,\
                     double p_relstop, double d_relstop, double rel_gap_relstop,\
                     double p_rho, double d_rho, double a_rho, double rho_g,\
                     double rhoI, double rhoM);

void calculate_residuals_no_structs(csi* Ai, csi* Aj, double* Av, csi m, csi n, csi nnz,
                                    double* b, double* c,\
                                    double* y, double* x, double tau, double* s, double kappa,\
                                    double p_relstop, double d_relstop, double g_relstop,\
                                    double* p_res, double* d_res, double* g_res,\
                                    double* n_p_res, double* n_d_res, double* n_g_res, double* rel_gap);
 

#ifdef __cplusplus
}
#endif

#endif

