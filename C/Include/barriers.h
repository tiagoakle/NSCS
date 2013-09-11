
#ifndef H_BARRIERS
#define H_BARRIERS
#include "common.h"
#include "nscs.h"

bool primal_feas(problem_t prob, double* x);
bool dual_feas(problem_t prob, double*x );
double barrier_complexity(problem_t prob);
void eval_grad(problem_t prob, double* x, double* grad);
void eval_hess(problem_t prob, double* x, state_t state);
csi hessian_nnz(problem_t prob);
csi cone_nnz(int type, csi n);
int cone_barrier_complexity(int type, csi n);
bool cone_dual_feas(int type, double* x, csi n);
bool cone_primal_feas(int type, double* x, csi n);
double cone_barrier_val(int type, double*x, csi n);
void cone_barrier_grad(int type, double*g, double*x, csi n);
csi cone_barrier_hessian(int type, int* HI, int* HJ, double* HV, double*x, csi n, double delta);
csi pos_orthant_complexity(csi n);
csi pos_orthant_nnz(int n);
double pos_orthant_val(double* x, csi n);
bool pos_orthant_feas(double* x, csi n);
void pos_orthant_grad(double* grad, double *x, csi n);
csi pos_orthant_hessian(int* HI, int* HJ, double* HV, double*x, csi n, double delta);


csi exp_complexity();
csi exp_nnz();
double exp_val(double* x);
bool exp_primal_feas(double* x);
bool exp_dual_feas(double* x);
void exp_grad(double* grad, double *x);
csi exp_hessian(int* HI, int* HJ, double* HV, double*x,double delta);

#endif
