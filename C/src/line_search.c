#include "nscs.h"
#include "nscs_sp.h"
#include "barriers.h"
#include "line_search.h"
#include "stdio.h"
#include "OpenBlAS/cblas.h"
#include "linear_solvers.h"
#include "eval_cent_meas.h"
#include <math.h>

/**
 * Backtracking linesearch for the aproximate tangent direction
 */
int linesearch_atd(state_t* state , parameters_t  pars, problem_t prob)
{
    double a0 = 1.; //Initial step length
    double a  = a0;

    // Calculate the largest step before either tau or kappa reach the boundary
    if(state->dkappa < 0) a = fmin(a,-state->kappa/state->dkappa);
    if(state->dtau < 0) a = fmin(a,-state->tau/state->dtau);
    //If with the full step either reaches the boundary make sure that after
    //the step the new trial tau or kappa is at most pars.eta of the way to the boundary.
    if(a<a0) a = a*pars.eta;

    //  allocate work vectors for the trial steps
    double* xa = (double*)calloc(state->x->n,sizeof(double));
    double* sa = (double*)calloc(state->x->n,sizeof(double));
    //Allocate work vectors to calculate the centrality
    double* psi  = (double*)calloc(state->x->n,sizeof(double));
    double* hpsi = (double*)calloc(state->x->n,sizeof(double));
    double kappaa;
    double taua;
    double dga;
    double mua;

    //TODO: add clean up before return?
    if(xa == NULL)   return OUT_OF_MEMORY;
    if(sa == NULL)   return OUT_OF_MEMORY;
    if(psi == NULL)  return OUT_OF_MEMORY;
    if(hpsi == NULL) return OUT_OF_MEMORY;
    
    //Counter for the number of backtracks
    int nsect = 0;
    int j = 0;

    //Define some variables
    //used in the loop.
    bool dFeas, pFeas; //Flags to indicate feasiblity
    bool dosect;       //Flag to indicate if there should be a backtrack

    double centmeas;   //Present value of centrality measure

    for(j=0;j<pars.max_backtrack;j++)
    {
          dosect = false;
    
          //Take an a sized step in x, tau and s, kappa, and store in xa, sa, taua, kappaa
          cblas_dcopy(state->x->n,state->x->v,1,xa,1);
          cblas_dcopy(state->s->n,state->s->v,1,sa,1);
          cblas_daxpy(prob.n,a,state->dx->v,1,xa,1);
          cblas_daxpy(prob.n,a,state->ds->v,1,sa,1);
          taua      = state->tau + a*state->dtau;
          kappaa    = state->kappa + a*state->dkappa;
            
          //Calculate the duality gap at the new trial point
          dga    = cblas_ddot(state->x->n,xa,1,sa,1) + taua*kappaa;
          mua    = dga / (prob.nu + 1); 
         
          //Check if x,s are feasible wrt the cones
          dFeas = dual_feas(prob,sa);
          pFeas = primal_feas(prob,xa);
    
         //If not primal or dual feasible backtrack
        if(!pFeas)
        {
            dosect = true; 
            //printf("C: Not primal feasible\n");
        }
        else if(!dFeas)
        {
            dosect = true; 
            //printf("C: Not dual feasible\n");
        }
        else //If both primal and dual feasible measure the centrality
        {
            centmeas = eval_cent_meas(prob,xa,sa,*state,mua,psi,hpsi);
            //Decide if we need to backtrack
            if(centmeas > mua*pars.theta)
            {
                dosect = true;
            }
            
            //printf("Evaluating dist: %g, %g, %i:\n",centmeas,mua*pars.theta,dosect);
        }
        //If we must backtrack do so
        if(dosect)
        {
            a = a*pars.lscaff;
            state->nbacktrack += 1; 
        }
        else
        {
            break;
        }
         
    }  
 
    //Check if the backtrack reached the maximum number of iterates
    if(j==pars.max_backtrack)
    {
       //If the backtrack fails do not update the vectors and release 
       //the local work vectors, and test vectors 
       free(xa);
       free(sa);
       free(psi);
       free(hpsi);
       state->a = a;   
       return BACKTRACK_FAIL; 
    }
    else
    {
        //If the backtrack succedded accept the iterate
        //store it in state and take the step in y

        //First release the old vectors
        double* xt = state->x->v;
        double* st = state->s->v;
        
        state->x->v    = xa;
        state->s->v    = sa;
        state->tau    = taua;
        state->kappa  = kappaa;
       
        //Take the step in y
        cblas_daxpy(prob.m,a,state->dy->v,1,state->y->v,1);
        state->a = a;   

        //Free the old vectors
        free(xt);
        free(st);
        free(psi);
        free(hpsi);

        return OK;
    }

    //state->a = a;   
    //return OK;
}


