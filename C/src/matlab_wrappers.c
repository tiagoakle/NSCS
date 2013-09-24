#include "linear_solvers.h"
#include "spmat.h"
#include "nscs.h"
#include "matlab_include/matlab_lib.h"
#include "eval_cent_meas.h"
//Contains the matlab wrappers for certain calls to the library

//The following functions should substitute the solve KKT system call
//function [d,CF] = solve_linear_system(H,mu,A,b,c,tau,kappa,r1,r2,r3,r4,r5,pars)

//Wraps solve_kkt_system so no structures are used and it is easier to call from 
//MATLAB

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
                                double* dk )
{
   //Build the structures
   spmat H;
   H.m   = n;
   H.n   = n;
   H.nnz = Hnnz;
   H.I   = Hi;
   H.J   = Hj;
   H.V   = Hv;
  
   spmat A;
   A.m   = m;
   A.n   = n;
   A.nnz = Annz;
   A.I   = Ai;
   A.J   = Aj;
   A.V   = Av;


double scalars[] = {m,n,tau,kappa};

//Now execute the call
int ret = solve_kkt_system(mu,\
                     H,\
                     A,\
                     b,\
                     c,\
                     tau,\
                     kappa,\
                     delta,\
                     gamma,\
                     r1,\
                     r2,\
                     r3,\
                     r4,\
                     r5,\
                     dy,\
                     dx,\
                     dt,\
                     ds,\
                     dk );

//Debug stuff
write_vector_to_csv("dx_call.csv", dx, n);

return ret;

}

/**
 * Wrapper to call from matlab
 */
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
                                                                                                    csi nnzH,\
                                                                                                    double* a,\
                                                                                                    int* nbacktrack)
{

    //Construct prob
    spmat A;
    A.n = n;
    A.m = m;
   
    problem_t prob;
    prob.A       = A;
    prob.tK      = tK;
    prob.nK      = nK;
    prob.k_count = k_count;
    prob.delta   = 1.e-10;
    prob.gamma   = 1.e-10;
    prob.nu      = nu;

    //Construct the state
    spmat H;
    H.nnz = nnzH;

    state_t state;
    state.H = H;
    state.y = y;
    state.x = x; 
    state.s = s;
    state.tau = tau;
    state.kappa = kappa;

    state.dy = dy;
    state.dx = dx;
    state.ds = ds;
    state.dtau = dtau;
    state.dkappa = dkappa;
    
    state.mu     = NAN; //In case we access this
    state.nbacktrack = 0;
    
    //Build the params structure

    parameters_t params;
    params.print = true;
    params.max_center_iter = NAN;
    params.eta   = eta;
    params.theta = theta;
    params.lscaff = lscaff;
    params.max_backtrack = max_backtrack;

    //Allocate space for the hessian!!!
    state.H.I   = (int*)calloc(nnzH,sizeof(int));
    state.H.J   = (int*)calloc(nnzH,sizeof(int));
    state.H.V   = (double*)calloc(nnzH,sizeof(double));
    state.H.n   = n;
    state.H.nnz = nnzH;

    //Call the linesearch
    linesearch_atd(&state,params,prob);
    (*nbacktrack) = state.nbacktrack;
    (*a)          = state.a;
    
    free(state.H.I);
    free(state.H.J);
    free(state.H.V);

}

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
                                                                                                    double* objvalb)
{

    //Construct prob
    spmat A;
    A.n = n;
    A.m = m;

    //Construct prob
    problem_t prob;
    prob.A       = A;
    prob.tK      = tK;
    prob.nK      = nK;
    prob.k_count = k_count;
    prob.delta   = 1.e-10;
    prob.gamma   = 1.e-10;
    prob.nu      = nu;

    //Construct the state
    spmat H;
    H.nnz = nnzH;

    state_t state;
    state.H = H;
    state.y = y;
    state.x = x; 
    state.s = s;
    state.tau = tau;
    state.kappa = kappa;
    state.mu    = mu;
    state.min_centmeas = *objvalb; //Before the call assign the value to the calling param value

    state.dy = dy;
    state.dx = dx;
    state.ds = ds;
    state.dtau = dtau;
    state.dkappa = dkappa;
    
    state.mu     = NAN; //In case we access this
    state.nbacktrack = 0;
    
    //Build the params structure

    parameters_t params;
    params.print = true;
    params.max_center_iter = NAN;
    params.eta   = eta;
    params.theta = theta;
    params.lsccent = lsccent;
    params.max_backtrack = max_backtrack;

    //Allocate space for the hessian!!!
    state.H.I   = (int*)calloc(nnzH,sizeof(int));
    state.H.J   = (int*)calloc(nnzH,sizeof(int));
    state.H.V   = (double*)calloc(nnzH,sizeof(double));
    state.H.n   = n;
    state.H.nnz = nnzH;

    //Call the linesearch
    linesearch_centering(&state,params,prob);
    //Return the results
    (*nbacktrack) = state.nbacktrack;
    (*a)          = state.a;
    (*objvalb)     = state.min_centmeas;
    
    free(state.H.I);
    free(state.H.J);
    free(state.H.V);

}

//Evaluates the hessian
void eval_hess_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* HI, int* HJ, double* HV)
{
   problem_t prob;
   prob.k_count = k_count;
   prob.tK = tK;
   prob.nK = nK;
   prob.delta = 0.0;
   state_t stat;
   spmat H;
   stat.H = H;
   stat.H.I = HI;
   stat.H.J = HJ;
   stat.H.V = HV;
   eval_hess(prob,x,stat);
}

//Evaluates the gradient
void eval_grad_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, double* g)
{
   problem_t prob;
   prob.k_count = k_count;
   prob.tK = tK;
   prob.nK = nK;
   prob.delta = 0.0;
   eval_grad(prob,x,g);
}

//Evaluates the hessian
void primal_feas_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* feas)
{
   problem_t prob;
   prob.k_count = k_count;
   prob.tK = tK;
   prob.nK = nK;
   prob.delta = 0.0;
   state_t stat;
   *feas = primal_feas(prob,x);
}

//Evaluates the dual feasibility
void dual_feas_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* feas)
{
   problem_t prob;
   prob.k_count = k_count;
   prob.tK = tK;
   prob.nK = nK;
   prob.delta = 0.0;
   state_t stat;
   *feas = dual_feas(prob,x);
}

void eval_cent_meas_no_structs(csi k_count, csi* nK, int* tK, double delta, double* x, double* s, csi n, csi nnzH, double mua, double* psi, double * hpsi, double* centmeas)
{
   //Construct prob
   spmat A;
   A.n = n;

   problem_t prob;
   prob.A       = A;
   prob.k_count = k_count;
   prob.tK = tK;
   prob.nK = nK;
   prob.delta = 0.0;
    //Allocate space for the hessian
   state_t state; 
   state.H.I   = (csi*)calloc(nnzH,sizeof(csi));
   state.H.J   = (csi*)calloc(nnzH,sizeof(csi));
   state.H.V   = (double*)calloc(nnzH,sizeof(double));
   state.H.n   = n;
   state.H.nnz = nnzH;
   //Call the function
   *centmeas = eval_cent_meas(prob,  x,  s, state,  mua, psi, hpsi);
   free(state.H.I);
   free(state.H.J);
   free(state.H.V);
}
