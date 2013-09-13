#include "linear_solvers.h"
#include "nscs_sp.h"
#include "nscs.h"
#include "matlab_include/matlab_lib.h"
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

   vec   vb;
   vb.n   = m;
   vb.v   = b;

   vec   vc;
   vc.n   = n;
   vc.v   = c;

   vec   vr1;
   vr1.n  = m;
   vr1.v  = r1;

   vec   vr2;
   vr2.n  = n;
   vr2.v  = r2;

   vec   vr4;
   vr4.n  = n;
   vr4.v  = r4;
   
   //Now build the vectors where the result will be stored
   //we assume that dy,dx,ds are already allocated
   vec   vdy;
   vdy.n  = m;
   vdy.v  = dy;

   vec   vdx;
   vdx.n  = n;
   vdx.v  = dx;

   vec    vds;
   vds.n  = n;
   vds.v  = ds;

double scalars[] = {m,n,tau,kappa};
//Debug stuff
//write_vector_to_csv("b_call.csv", vb.v, vb.n);
//write_vector_to_csv("c_call.csv", vc.v, vc.n);
//write_vector_to_csv("r1_call.csv", vr1.v, vr1.n);
//write_vector_to_csv("r2_call.csv", vr2.v, vr2.n);
//write_vector_to_csv("r4_call.csv", vr4.v, vr4.n);
//write_vector_to_csv("doubles_call.csv",scalars,4);

//Now execute the call
int ret = solve_kkt_system(mu,\
                     H,\
                     A,\
                     vb,\
                     vc,\
                     tau,\
                     kappa,\
                     delta,\
                     gamma,\
                     vr1,\
                     vr2,\
                     r3,\
                     vr4,\
                     r5,\
                     vdy,\
                     vdx,\
                     dt,\
                     vds,\
                     dk );

//Debug stuff
write_vector_to_csv("dx_call.csv", vdx.v, vdx.n);

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
    problem_t prob;
    prob.tK      = tK;
    prob.nK      = nK;
    prob.k_count = k_count;
    prob.m       = m;
    prob.n       = n;
    prob.delta   = 1.e-10;
    prob.gamma   = 1.e-10;
    prob.free    = 0;
    prob.nu      = nu;
    prob.nnzH    = nnzH;

    //Construct the state
    state_t state;
    vec vy;
    state.y = &vy;
    state.y->v = y;
    state.y->n = m;

    vec vx;
    state.x = &vx;
    state.x->v = x;
    state.x->n = n;
    
    vec vs;
    state.s = &vs;
    state.s->v = s;
    state.s->n = n;
    
    state.tau = tau;
    state.kappa = kappa;


    vec vdy;
    state.dy = &vdy;
    state.dy->v = dy;
    state.dy->n = m;

    vec vdx;
    state.dx = &vdx;
    state.dx->v = dx;
    state.dx->n = n;
    
    vec vds;
    state.ds = &vds;
    state.ds->v = ds;
    state.ds->n = n;
    
    state.dtau = dtau;
    state.dkappa = dkappa;
    
    state.mu     = NAN; //In case we access this
    state.eta    = NAN;
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
    state.H.I = (int*)calloc(nnzH,sizeof(int));
    state.H.J = (int*)calloc(nnzH,sizeof(int));
    state.H.V = (double*)calloc(nnzH,sizeof(double));

    //Call the linesearch
    linesearch_atd(state,params,prob);
    (*nbacktrack) = state.nbacktrack;
    (*a)          = state.a;
    
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
