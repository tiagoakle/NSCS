#ifndef H_LINE_SEARCH
#define H_LINE_SEARCH
//int linesearch_atd(state_t state ,parameters_t  pars, problem_t prob);
#ifdef __cplusplus
extern "C" {
#endif

int linesearch_atd_no_structs( int m, int n, double*x, double*y, double*s, double tau, double kappa, double* dx,\
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
#ifdef __cplusplus
}
#endif


#endif
