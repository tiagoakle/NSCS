#include "eval_cent_meas.h"
#include "nscs.h"
#include "stdio.h"
#include "OpenBLAS/cblas.h"
//Evaluates the centrality measure, returns grad to debug only
double eval_cent_meas(problem_t prob, double* xa, double* sa, state_t state, double mua, double* psi, double * hpsi)
{

        //Evaluate the gradient at the preset point
        eval_grad(prob,xa,psi); 
        eval_hess(prob,xa,state);
        //scale by mu and add s to psi
        cblas_dscal(prob.n,mua,psi,1);
        cblas_daxpy(prob.n,1.0,sa,1,psi,1);
        //Solve the linear system H(hpsi)=psi
        int ret = solve_linear_system(hpsi,state.H.I,state.H.J,state.H.V,state.H.nnz,psi,state.H.n);
        if(ret!=0){printf("Error solving linear system ret:%i, n:%i, nnz:%i",ret,state.H.n, state.H.nnz); return INTERNAL_ERROR; }//XXX:We should check for numerical error and not crash 
        double centmeas = sqrt(cblas_ddot(prob.n,psi,1,hpsi,1)); 
        return centmeas;
}
