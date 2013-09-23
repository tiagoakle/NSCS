#include <math.h>
#include <stdio.h>
#include "nscs.h"
#include "spmat.h"
#include "barriers.h"
#include "line_search.h"
#include "OpenBlAS/cblas.h"
#include "linear_solvers.h"
#include "eval_cent_meas.h"

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
    double* xa = (double*)calloc(prob.A.n,sizeof(double));
    double* sa = (double*)calloc(prob.A.n,sizeof(double));
    //Allocate work vectors to calculate the centrality
    double* psi  = (double*)calloc(prob.A.n,sizeof(double));
    double* hpsi = (double*)calloc(prob.A.n,sizeof(double));
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
          cblas_dcopy(prob.A.n,state->x,1,xa,1);
          cblas_dcopy(prob.A.n,state->s,1,sa,1);
          cblas_daxpy(prob.A.n,a,state->dx,1,xa,1);
          cblas_daxpy(prob.A.n,a,state->ds,1,sa,1);
          taua      = state->tau + a*state->dtau;
          kappaa    = state->kappa + a*state->dkappa;
            
          //Calculate the duality gap at the new trial point
          dga    = cblas_ddot(prob.A.n,xa,1,sa,1) + taua*kappaa;
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

        //Copy the successeful xa and sa to state.x state.s
        //just switching the pointers does not work because
        //we should not free a matlab allocated vector
        cblas_dcopy(prob.A.n,xa,1,state->x,1);
        cblas_dcopy(prob.A.n,sa,1,state->s,1);
        state->tau    = taua;
        state->kappa  = kappaa;
       
        //Take the step in y
        cblas_daxpy(prob.A.m,a,state->dy,1,state->y,1);
        state->a = a;   

        //Free the work vectors
        free(xa);
        free(sa);
        free(psi);
        free(hpsi);

        return OK;
    }

    //state->a = a;   
    //return OK;
}

/**
 * Backtracking linesearch for the centering direction
 */
int linesearch_centering(state_t* state , parameters_t  pars, problem_t prob)
{
    double a0 = 1.; //Initial step length
    double a  = a0*pars.eta;

    // Calculate the largest step before either tau or kappa reach the boundary
    if(state->dkappa < 0) a = fmin(a,-state->kappa/state->dkappa);
    if(state->dtau < 0) a = fmin(a,-state->tau/state->dtau);
    //If with the full step either reaches the boundary make sure that after
    //the step the new trial tau or kappa is at most pars.eta of the way to the boundary.
    if(a<a0) a = a*pars.eta;

    //  allocate work vectors for the trial steps
    double* xa = (double*)calloc(prob.A.n,sizeof(double));
    double* sa = (double*)calloc(prob.A.n,sizeof(double));
    //Allocate work vectors to calculate the centrality
    double* psi  = (double*)calloc(prob.A.n,sizeof(double));
    double* hpsi = (double*)calloc(prob.A.n,sizeof(double));
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
    double min_centmeas = state->min_centmeas; //Minimum centrality found 

    for(j=0;j<pars.max_backtrack;j++)
    {

          //printf("C: backtrack: nsect %i :a %g, objval: %g \n",nsect,a,min_centmeas);
          dosect = false;
    
          //Take an a sized step in x, tau and s, kappa, and store in xa, sa, taua, kappaa
          cblas_dcopy(prob.A.n,state->x,1,xa,1);
          cblas_dcopy(prob.A.n,state->s,1,sa,1);
          cblas_daxpy(prob.A.n,a,state->dx,1,xa,1);
          cblas_daxpy(prob.A.n,a,state->ds,1,sa,1);
          taua      = state->tau + a*state->dtau;
          kappaa    = state->kappa + a*state->dkappa;
            
          //Check if x,s are feasible wrt the primal cone only,
          //TODO: Check this statement
          //the feasiblity of the dual is assured by the objective 
          //centrality 
          pFeas = primal_feas(prob,xa);
    
         //If not primal feasible backtrack
        if(!pFeas)
        {
            dosect = true; 
            //printf("C: Not primal feasible\n");
        }
        else //If primal feasible measure the centrality
        {
            centmeas = eval_cent_meas(prob,xa,sa,*state,state->mu,psi,hpsi);
            //If the centrality here is smaller than the minimum yet then backtrack
            //if it is not, then the centrality measure started increasing 
            if(centmeas < min_centmeas)
            {
                min_centmeas = centmeas;
                dosect = true;
            }
            
            //printf("Evaluating dist: %g, %g, %i:\n",centmeas,mua*pars.theta,dosect);
        }
        
        //If we must backtrack do so
        if(dosect)
        {
            a = a*pars.lsccent;
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

        //Copy the successeful xa and sa to state.x state.s
        //just switching the pointers does not work because
        //we should not free a matlab allocated vector
        cblas_dcopy(prob.A.n,xa,1,state->x,1);
        cblas_dcopy(prob.A.n,sa,1,state->s,1);
        state->tau    = taua;
        state->kappa  = kappaa;
       
        //Take the step in y
        cblas_daxpy(prob.A.m,a,state->dy,1,state->y,1);
        state->a = a;   
        state->min_centmeas = min_centmeas;//Store the centrality in the state structure

        //Free the work vectors
        free(xa);
        free(sa);
        free(psi);
        free(hpsi);

        return OK;
    }

    //state->a = a;   
    //return OK;
}
