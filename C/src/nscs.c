#include "nscs.h"
#include "spmat.h"
#include "barriers.h"
#include "common.h"
#include "smatvec.h"
#include "OpenBLAS/cblas.h"
#include "line_search.h"
#include "linear_solvers.h"
#include <stdio.h>
#include <math.h>

/**
 * Main entry point for the Non symmetric cone solver.
 *
 *
 * NSCS solves a conic programming problem given in 
 * standard form,
 * min c'x st Ax=b x \in K.
 * K is the product of several cones defined by the
 * arrays problem.nK, problem.iK. 
 * 
 * K = K_0xK_1xK_2,.....,K_q-1
 *
 * Where 
 * each cone is one of: positive orthant (0), second order cone (1), 
 * semi-definite cone (2), exponential cone (3) , power cone (4).
 * 
 * tK is an array of integers and of length q. The value tK[0] indicates the type
 * of cone 0,..., tK[q-1] indicates the type of cone K_{q-1}.
 *
 * nK is a length q array of integers. It parametrizes the size of the cones,
 * nK[0] the size of cone 1,..., nK[k_cones-1] the size 
 * of cone k_cones-1. All exponential cones and power cones must be of size 3. Their
 * corresponding nK[.] cannot be ignored and must be set to 3.
 *
 * @param problem contains the problem definition 
 * @param pars    contains the algorithm parameters
 * @param result  is populated with the results
 * @param x0      initial primal point, mandatory
 * @param y0      initial dual y not mandatory, can be set to null
 * @param s0      initial dual s not mandatory, can be set to null
 * @param tau     initial tau value, not mandatory can be set to null
 * @param kappa   initial kappa value, not mandatory can be set to null
 * @return status flag
 *
 */
int nscs(problem_t* problem, parameters_t* pars,\
        double* y, double* x, double* t, double* s, double* k, bool wy, bool wt, bool ws, bool wk) 
{
    
    print_header();
    
    //Define some suporting variables
    int status;
    
    //Validate the parameters
    status = validate_pars(problem,pars);
    if(status != OK)
    {
        fprintf(stderr,"invalid parameter provided");
        return status;

    } 

    //Initalize the state structure
    state_t state;    
    //Allocate the structure to hold the state of the algorithm
    status     =  allocate_state(&state, problem);
    if(status != OK)
    {
        free_state(state);
        return status;
    }

    //Initialize the variables that are a function of the parameters and
    //problem definition
    status = init(problem,&state,pars,y,x,t,s,k,wy,wt,ws,wk);
    if(status != OK)
    {
        fprintf(stderr,"error initializing algorithm");
        return status;
    }

    //Initialize the residuals
    calculate_residuals(&state,problem);
   
    //Flag that indicates that extra centering steps are needed
    bool do_center;
    int c_iter;
    //Integer to hold the reason to stop 
    int stop_reason = END_CONTINUE;
    
    //Start of main loop ---------------
    for(state.m_iter = 0;state.m_iter < pars->max_iter; state.m_iter++)
    {
        //Evaluate the barrier function
        eval_hess(*problem,state.x,state);
        //Calculate the approximate direction
        status = solve_approximate_tangent_direction(&state,*problem,*pars);
        
        //Save the previous tau and kappa
        state.prev_tau = state.tau;
        state.prev_kappa = state.kappa;
        //Do a line search 
        status = linesearch_atd(&state,*pars,*problem);
        if(status!=OK)
        {
            fprintf(stderr,"Approximate tangent direction linesearch failed");
            free_state(state);
            return status;
        }
        //Calculate the new residuals
        calculate_residuals(&state,problem);
        //Update the centrality measurements
        state.dgap  = cblas_ddot(problem->A.n,state.x,1,state.s,1) + state.tau*state.kappa;
        state.mu    = state.dgap/(problem->nu+1);
        state.vfeas = state.dtau/state.prev_tau - state.dkappa/state.prev_kappa;
      

        //At this point the approximate tangent direction linesearch is finished now 
        //start the centering direction
        c_iter = 0;
        //Run the centering routine
        while(do_center && c_iter < pars->max_center_iter)
        { 
            //Calculate the centering direction and set the value
            //of dx_norm
            status = calculate_centering_direction(&state,*problem, pars);
            //Linesearch for centering 
            status = linesearch_centering(&state , *pars, *problem);
            if(status!=OK)
            {
               fprintf(stderr,"Approximate tangent direction linesearch failed");
               free_state(state);
               return status;
            }
            c_iter++;

            //Check if the centering condition is satisfied
            if(state.dx_norm<pars->beta)
            {
                do_center = false;
            }
        }
        state.ncent_iter = c_iter;

        //If the maximum number of iterations was exceeded then exit with an error
        if(do_center)
        {
            stop_reason = END_MAX_CENTER_ITER;
            break;
        }

        //Print the output and check for exit
        if(pars->print)
        { 
            print_and_log(&state,pars);
        }
       
        //Check the stopping criteria
        stop_reason = check_stopping_criteria(&state,pars,problem);
        //If stop_reason is not END_CONTINUE, then
        //the iteration must stop
        if(stop_reason!=END_CONTINUE)
        {
            break;
        }

    } //End of the main loop-----------------------------
    
    if(state.m_iter == pars->max_iter)
    {
        stop_reason = END_MAX_ITER;
    }
    
    //If optimal scale the results
    if(stop_reason == END_OPTIMAL)
    {
        cblas_dscal(problem->A.n,1./state.tau,state.x,1);
        cblas_dscal(problem->A.n,1./state.tau,state.s,1);
        cblas_dscal(problem->A.m,1./state.tau,state.y,1);
    }

    //Print final info
    if(pars->print)
    {
        print_final(stop_reason);
    }

    //Clean up
    free_state(state);
    return 0;
}

/**
 * Function to free the state structure
 * It will not free y,x,s 
 * so as to keep them for the result strucutre. 
 */
void free_state(state_t state)
{
    //Free the present directions
    if(!(state.dy==NULL))   free(state.dy);
    if(!(state.dx==NULL))   free(state.dx);
    if(!(state.ds==NULL))   free(state.ds);
   
    //(Free the residual vectors
    if(!(state.p_res==NULL))   free(state.p_res);
    if(!(state.d_res==NULL))   free(state.d_res);
    //if(!(state.c_res==NULL))   free(state.c_res);
   
}

/**
 * Validates the parameters provided by the user
 */
int validate_pars(problem_t* problem, parameters_t* pars)
{
    if(val_range(pars->max_iter,1,1000,"max_iter")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->max_center_iter,1,1000,"max_center_iter")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->theta,1,1000,"theta")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->lscaff,1.e-2,1.-1.e-2,"lscaff")!=OK) return INVALID_PARAMETER;;
    if(val_range(pars->lsccent,1.e-2,1.-1.e-2,"lsccent")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->eta,0.5,1.-1.e-10,"eta")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->beta,0.1,1.-1.e-10,"beta")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->p_relstop,1.e-16,1.e-5,"p_relstop")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->d_relstop,1.e-16,1.e-5,"d_relstop")!=OK) return INVALID_PARAMETER;
    if(val_range(pars->rel_gap_relstop,1.e-16,1.e-5,"rel_gap_relstop")!=OK) return INVALID_PARAMETER;
    //Check the cones
    //Check that the cone is defined 
    return VALIDATION_OK;
}

int val_range(double param, double min, double max, char* name)
{
    if(param<min||param>max)
    { 
        fprintf(stderr,"parameter %s: %g outside range %g, %g",name,param,min,max);
        return INVALID_PARAMETER;
    }
    return OK;
}

//Calculate the starting iterate and save it in state_t
int allocate_state(state_t* state,problem_t* problem)
{
    
    //Allocate the vectors that hold the present residual
    state->p_res = (double*)calloc(problem->A.m,sizeof(double));
    if(state->p_res == NULL) return OUT_OF_MEMORY;
    state->d_res = (double*)calloc(problem->A.n,sizeof(double));
    if(state->d_res == NULL) return OUT_OF_MEMORY;
    //state->c_res = (double*)calloc(problem->A.n,sizeof(double));
    //if(state->c_res == NULL) return OUT_OF_MEMORY;

    //Allocate the vectors that hold the present search direction
    state->dy = (double*)calloc(problem->A.m,sizeof(double));
    if(state->dy == NULL) return OUT_OF_MEMORY;
    state->dx = (double*)calloc(problem->A.n,sizeof(double));
    if(state->dx == NULL) return OUT_OF_MEMORY;
    state->ds = (double*)calloc(problem->A.n,sizeof(double));
    if(state->ds == NULL) return OUT_OF_MEMORY;
    
    //allocate space for the stateient
    state->grad = (double*)calloc(problem->A.n,sizeof(double));
    if(state->grad == NULL) return OUT_OF_MEMORY;
    
    return OK;
}

/**
 * Calculates the residuals and sets the appropriate 
 * variables in state.
 * @param state pointer to the state structure
 * @param prob  problem definition stucture
 */
void calculate_residuals(state_t* state, problem_t* prob)
{
     //-Ax+tb
     cblas_dcopy(prob->A.m,prob->b,1,state->p_res,1);
     cblas_dscal(prob->A.m,state->tau,state->p_res,1); 
     //Product in coordinate form
     dspmvcoo(prob->A.nnz, -1.0, prob->A.I, prob->A.J,prob->A.V, state->p_res, state->x);

     //A'y+s-tc: 
     cblas_dcopy(prob->A.n,state->s,1,state->d_res,1);
     dsptmvcoo(prob->A.nnz,1., prob->A.I, prob->A.J, prob->A.V, state->d_res, state->y);
     cblas_daxpy(prob->A.n,-state->tau,prob->c,1,state->d_res,1);
     //Transpose product in coordinate form

     //-b'y+c'x+k
     double bty = cblas_ddot(prob->A.m,prob->b,1,state->y,1);
     double ctx = cblas_ddot(prob->A.n,prob->c,1,state->x,1);
       state->g_res = -bty
                     +ctx\
                     +state->kappa; 
   
    //Relative_gap
    state->rel_gap = fabs(ctx-bty)/fabs(state->tau+fabs(bty)); 
    //infinity norm of p_res
    csi ix = cblas_idamax(prob->A.m,state->p_res,1);
    state->n_p_res = fabs(state->p_res[ix])/state->p_relstop;
    //Relative infinity norm of d_res
    ix = cblas_idamax(prob->A.n,state->d_res,1);
    state->n_p_res = fabs(state->d_res[ix])/state->d_relstop;
    //Relative infinity norm of 
    state->n_g_res   = state->g_res/state->g_relstop;

}


/**
 * This solves for the approximate tangent direction 
 * and stores it in state
 */
int  solve_approximate_tangent_direction(state_t *state,problem_t prob,parameters_t pars)
{
//[d,CF] = solve_linear_system(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,-v.rP,-v.rD,-v.rG,-v.kappa*v.tau,-v.s,pars); 
  //TODO: Change the definition of the residuals so that they match the parameter of the function
  //Define the rhs 
  //Allocate space for -s
  double* ms = (double*)calloc(prob.A.n,sizeof(double));
  cblas_daxpy(prob.A.n,-1.,state->s,1,ms,1); //ms = -s

  solve_kkt_system(state->mu,\
                   state->H,\
                   prob.A,\
                   prob.b,\
                   prob.c,\
                   state->tau,\
                   state->kappa,\
                   pars.delta,\
                   pars.gamma,\
                   state->p_res,\
                   state->d_res,\
                   state->g_res,\
                   ms,\
                   -state->kappa*state->tau,\
                   state->dy,\
                   state->dx,\
                   &state->dtau,\
                   state->ds,\
                   &state->dkappa);
    free(ms);
    return OK;
}

/**
 * Calculates the rhs for the linear system that defines
 * the centering direction and evaluates ||dx||_x < beta.
 * If this is true it sets do_center to false.
 */
int calculate_centering_direction(state_t *state, problem_t prob, parameters_t* pars)
{
    //            % right-hand sides:
    //        rs.r1 = zeros(size(pars.A,1),1);
    //        rs.r2 = zeros(size(pars.A,2),1);
    //        rs.r3 = 0;
    //        rs.r4 = v.mu   - (v.tau*v.kappa);
    //        rs.r5 = - v.s-v.mu*v.F{2}; %        - (-pars.A'*v.y + v.mu*v.F{2} + v.tau*pars.c);
    //Allocate space for the zeros
    double* rs1 = (double*)calloc(prob.A.m,sizeof(double));
    double* rs2 = (double*)calloc(prob.A.n,sizeof(double));
    double* cent = (double*)calloc(prob.A.n,sizeof(double));
    eval_grad(prob,state->x,cent);
    cblas_dscal(prob.A.n,-state->mu,cent,1);
    cblas_daxpy(prob.A.n,-1.,state->s,1,cent,1);

    int res =  solve_kkt_system(state->mu,\
                   state->H,\
                   prob.A,\
                   prob.b,\
                   prob.c,\
                   state->tau,\
                   state->kappa,\
                   pars->delta,\
                   pars->gamma,\
                   rs1,\
                   rs2,\
                   0.,\
                   cent,\
                   state->mu -state->kappa*state->tau,\
                   state->dy,\
                   state->dx,\
                   &state->dtau,\
                   state->ds,\
                   &state->dkappa);

    //Use cent, to evaluate ||dx||_x
    cblas_daxpy(prob.A.n,-1.,state->ds,1,cent,1);
    state->dx_norm = sqrt(cblas_ddot(prob.A.n,state->dx,1,cent,1)/state->mu);
    free(rs1);
    free(rs2);
    free(cent); 
}

/**
 * Print the present state of the iteration
 *
 */
void print_and_log(state_t* state,parameters_t* pars)
{
    fprintf(stdout,"%i  %i  %e  %e  %e  %e  %e  %e  %e  %e\n",state->m_iter,state->ncent_iter,state->a,state->mu,state->tau\
                                                             ,state->kappa,state->n_p_res,
                                                             state->n_d_res,state->rel_gap,state->n_g_res);
}
/**
 * Print header
 *
 */
void print_header()
{
    fprintf(stdout,"======================================\n");
    fprintf(stdout,"%s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n",
    "iter","center","a","mu","tau",\
    "kappa","p_res","d_res","rel_gap","g_res");
    fprintf(stdout,"======================================\n");
}

void print_final(int stop_reason)
{
    char* reason;
    switch(stop_reason)
    {
        case END_OPTIMAL:
            reason = "Optimal";
        break;
        case END_D_INFEAS:
            reason = "Dual infeasible";
        break;
        case END_P_INFEAS:
            reason = "Primal infeasible";
        break;
        case END_ILL_POSED:
            reason = "Ill Posed";
        break;
        case END_MAX_ITER:
            reason = "Max iter reached";
        break;
        case END_MAX_CENTER_ITER:
            reason = "Max number of center iterations reached"; 
        break;
    }
    fprintf(stdout,"======================================\n");
}

//Check stopping criteria
int check_stopping_criteria(state_t *state,parameters_t* pars, problem_t* prob)
{
    //Check if all the residuals are smaller than the tolearances
    if(state->n_p_res<pars->p_relstop)
        if(state->n_d_res<pars->d_relstop)
            if(state->rel_gap<pars->rel_gap_relstop) //Optimal solution
            {
               return END_OPTIMAL; 
            }
            else
            {
                if(state->n_g_res<pars->rho_g) //check if ctx-bty-k like 0
                {
                    if(state->tau<pars->rhoI*fmax(1.,state->kappa)) //Small relative gap and tau small wrt kappa infeasible
                    {
                        bool p_infeas, d_infeas;
                        double ctx = cblas_ddot(prob->A.n,state->x,1,prob->c,1);
                        double bty = cblas_ddot(prob->A.m,state->y,1,prob->b,1);
                        d_infeas = (ctx<0);
                        p_infeas = (bty>0);
                        
                        if(p_infeas&&d_infeas)//Do another comparisson and decide
                        {
                           if((fabs(ctx)-fabs(bty))>0) //ctx has larger magnitude
                           { 
                                d_infeas = true;
                                p_infeas = false; 
                                return END_D_INFEAS;
                           }
                           else //bty has larger magnitude
                           {
                                d_infeas = false;
                                p_infeas = true; 
                                return END_P_INFEAS;
                           }
                        }
                        else if(p_infeas)
                        {
                                return END_P_INFEAS;
                        }
                        else if(d_infeas)
                        {
                                return END_D_INFEAS;
                        }
                    }
                }
            }

   //XXX: This ill-posedness check seems wrong 
   if(state->mu<pars->rhoM*state->mu0)
        if(state->tau<pars->rhoI*fmin(1,state->kappa))
        {
            return END_ILL_POSED;
        }
}

/**
 * Initializes the variables that are a function of the parameters and the initial points
 * except for the residuals.
 * Counts the nnz of H and allocates space for state.H, sets state.H.nnz state.H.n state.H.m
 * Sets the complexity parameter prob->nu
 * Calculates the initial point, if y,s,t,k are not provided
 * Calculates the initial value of mu, sets state->mu
 *
 * @param prob, the problem structure
 * @param state, the state structure
 * @param y0,t0,s0,k0 optional initial varaibles
 * @param x0          mandatory inital variables
 */
int init(problem_t* prob, state_t *state, parameters_t *pars, double* y0, double* x0, double* t0, double* s0, double* k0,\
         bool wy, bool wt, bool ws, bool wk)
{
     //Count the number of non zeros in the new hessian
    int k = 0;
    csi nnzH = 0; 
    for(k=0;k<prob->k_count;k++)
    {
        nnzH += cone_nnz(prob->tK[k], prob->nK[k]);
    }
    //Allocate the new hessian
    int status = calloc_spmat(&(state->H),prob->A.n,prob->A.n,nnzH);    
    if(status != OK) return status;
    //Initialize the complexity paramter
    prob->nu = barrier_complexity(*prob);

    //Define the initial point,
    //all the fields must be initialized, since they will be used as working vectors
    if(x0==NULL){ fprintf(stderr,"Initial point not provided \n"); return INVALID_PARAMETER;}
    if(y0==NULL){ fprintf(stderr,"Initial y has to be allocated \n"); return INVALID_PARAMETER;}
    if(s0==NULL){ fprintf(stderr,"Initial s has to be allocated \n"); return INVALID_PARAMETER;}
    if(t0==NULL){ fprintf(stderr,"Initial tau has to be allocated \n"); return INVALID_PARAMETER;}
    if(k0==NULL){ fprintf(stderr,"Initial kappa has to be allocated \n"); return INVALID_PARAMETER;}

    //Check if the primal is feasible
    if(!primal_feas(*prob,x0)){ fprintf(stderr,"Initial point not primal feasible\n"); return INVALID_PARAMETER;}
    //Check the dual is feasible if provided
    if(ws&&!dual_feas(*prob,s0)){ fprintf(stderr,"Initial dual not feasible\n"); return INVALID_PARAMETER;}
    //Check tau and kappa if provided
    if(wk&&*k0>0){ fprintf(stderr,"Initial kappa must be positive\n"); return INVALID_PARAMETER;}
    if(wt&&*t0>0){ fprintf(stderr,"Initial tau   must be positive\n"); return INVALID_PARAMETER;} 
    //complete the inital point
    state->x = x0;
    state->s = s0;
    state->y = y0;

    //Make sure tau and kappa are set
    if(wk)
    {
        state->kappa = *k0;
    }
    else
    {
        state->kappa = 1.0;
    }
    if(wt)
    {
        state->tau= *t0;
    }
    else
    { 
        state->tau = 1./state->kappa; //XXX:This is what coneopt does but does it make sense?
    }
    
    //Evaluate the initial gradient
    eval_grad(*prob,state->x,state->grad);
    if(!ws)
    {
        
        //As in coneopt
        // %v.s = -v.mu*v.F{2};
        
        // choosing s so that
        // s + (x'*s+tau*kappa)/(nu+1) * grad F(x) = 0
        // this is just a set of linear equations
        // solve them for s -- everything else is known
        // the system is easy to solve because
        // the system matrix is (I+v*x'). Use
        // Woodbury formula:
        // ---------------------------------------
        //vt  = v.F{2}/(K.nu+1);
        //qt  = -v.tau*v.kappa*v.F{2}/(K.nu+1);
        //v.s = qt-(v.x'*qt)*vt/(1+v.x'*vt);
        double* vt = (double*)calloc(prob->A.n,sizeof(double));
        if(vt==NULL){ fprintf(stderr,"Out of memory\n"); return OUT_OF_MEMORY;};

        cblas_daxpy(prob->A.n,1./(prob->nu+1.),state->grad,1,vt,1);
        //Calculate qt in s
        cblas_dcopy(prob->A.n,state->grad,1,state->s,1);
        cblas_dscal(prob->A.n,-state->tau*state->kappa/(prob->nu+1.),state->s,1);

        //Calculate the dot products
        double xqt = cblas_ddot(prob->A.n,state->x,1,state->s,1);
        double xvt = cblas_ddot(prob->A.n,state->x,1,vt,1);

        //Calculate the final s
        cblas_daxpy(prob->A.n,-xqt/(1+xvt),vt,1,state->s,1); 
        free(vt);

       
    }

    if(!wy)
    {
        int i;
        for(i=0;i<prob->A.m;i++)
        {
           state->y[i] = 1.0; 
        }
    }

    //Initialize mu and dgap
    state->dgap  = cblas_ddot(prob->A.n,state->x,1,state->s,1) + state->tau*state->kappa;
    state->mu    = state->dgap/(prob->nu+1);
    state->mu0   = state->mu;
    state->vfeas = 0.;

    //Initialize tauprev and kappa prev
    state->prev_tau = 1.;
    state->prev_kappa = 1.;
   
    //Initialize the counters
    state->ncent_iter = 0;

    //Copy the regularization parameters from parameters to problem
    prob->delta = pars->delta;
    prob->gamma = pars->gamma;
    
    return OK;
}

