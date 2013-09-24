#include "nscs.h"
#include "spmat.h"
#include "barriers.h"
#include "common.h"
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
int nscs(problem_t* problem, parameters_t* pars, result_t* result, double* y0, double* x0, double* t0, double* s0, double* k0) 
{

    //Define some suporting variables
    int status;
    
    //Validate the parameters
    status = validate_pars(problem,pars);
    if(status != OK)
    {
        fprintf(stderr,"invalid parameter provided");
        return status;

    } 

    //This variable holds a reference to the state structure
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
    status = init(problem,&state,y0,x0,t0,s0,k0);
    if(status != OK)
    {
        fprintf(stderr,"error initializing algorithm");
        return status;
    }

    //Main loop --------------------------------------------
    int m_iter = 0; //Main iteration counter
    
    bool do_center;
    int c_iter;
    for(m_iter = 0;m_iter < pars->max_iter; m_iter++)
    {
        //Evaluate the barrier function
        eval_hess(*problem,state.x,state);
        //Calculate the approximate direction
        status = solve_approximate_tangent_direction(state);
        //Do a line search 
        status = linesearch_atd(state,*pars,*problem);
        //Calculate the new residuals
        status = calculate_residuals(state);
        //Check the centering 
        do_center = !check_centering_condition(state,pars);
        c_iter = 0;
        //Run the centering routine
        while(do_center && c_iter < pars->max_center_iter)
        {
            status = calculate_centering_direction(state,pars);
            //Linesearch for centering 
            status = centering_linesearch(state,pars);
            //Check if the centering is complete
            do_center = !check_centering_condition(state,pars);
            c_iter++;
        }

        //If the maximum number of iterations was exceeded then clean and exit with 
        //an error
        if(do_center)
        {
            build_result(state,result,pars);
            free_state(state);
            return CENTER_ITERS_EXCEEDED;
        }

        //Print the output and check for exit
        if(pars->print)
        { 
            print_and_log(state,m_iter,c_iter,pars);
        }
       
       //Check the stopping criteria
       status = check_stopping_criteria(state,pars);

    } //End of the main loop-----------------------------

    //Print final info
    if(pars->print)
    {
        print_final(state);
    }

    build_result(state,result,pars);
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
    if(!(state.c_res==NULL))   free(state.c_res);
   
    if(!((state.H.I==NULL)&&(state.H.J==NULL)&&(state.H.V==NULL))) free_spmat(state.H);

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
    //Check the cones
    //Check that the cone is defined 
    return VALIDATION_OK;
}

int val_range(double param, double min, double max, char* name)
{
    if(param<min||param>max){ fprintf(stderr,"parameter %s: %g outside range %g, %g",name,param,min,max); return INVALID_PARAMETER;}
    return OK;
}

//Calculate the starting iterate and save it in state_t
int allocate_state(state_t* state,problem_t* problem)
{
    
    //Allocate the vectors that hold the present iterate
    state->y = (double*)calloc(problem->A.m,sizeof(double));
    if(state->y == NULL) return OUT_OF_MEMORY;
    state->x = (double*)calloc(problem->A.n,sizeof(double));
    if(state->x == NULL) return OUT_OF_MEMORY;
    state->s = (double*)calloc(problem->A.n,sizeof(double)); //s is the size of the dimensions of all constrained variables
    if(state->s == NULL) return OUT_OF_MEMORY;

    //Allocate the vectors that hold the present residual
    state->p_res = (double*)calloc(problem->A.m,sizeof(double));
    if(state->p_res == NULL) return OUT_OF_MEMORY;
    state->d_res = (double*)calloc(problem->A.n,sizeof(double));
    if(state->d_res == NULL) return OUT_OF_MEMORY;
    state->c_res = (double*)calloc(problem->A.n,sizeof(double));
    if(state->c_res == NULL) return OUT_OF_MEMORY;

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


int calculate_residuals(state_t state)
{
    return OK;
}

/**
 * This solves for the approximate tangent direction 
 * and stores it in state
 */
int  solve_approximate_tangent_direction(state_t state)
{
    return OK;
}

/**
 *
 * @return false if the iterate is not well centered
 */  
bool check_centering_condition(state_t state, parameters_t* pars)
{
    return true;
}

/**
 *
 */
int calculate_centering_direction(state_t state, parameters_t* pars)
{

}

int centering_linesearch(state_t state, parameters_t* pars)
{
}

void print_and_log(state_t state ,int m_iter, int c_iter,parameters_t* pars)
{
}

int check_stopping_criteria(state_t state,parameters_t* pars)
{
}

void print_final(state_t state)
{
}

void build_result(state_t state ,result_t* res, parameters_t* pars)
{
}


int init(problem_t* prob, state_t *state, double* y0, double* x0, double* t0, double* s0, double* k0)
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
    //Define the initial point
    if(x0==NULL){ fprintf(stderr,"Initial point not provided \n"); return INVALID_PARAMETER;}
    //Check if the primal is feasible
    if(!primal_feas(*prob,x0)){ fprintf(stderr,"Initial point not primal feasible\n"); return INVALID_PARAMETER;}
    //Check the dual is feasible if provided
    if(s0!=NULL&&!dual_feas(*prob,s0)){ fprintf(stderr,"Initial dual not feasible\n"); return INVALID_PARAMETER;}
    //Check tau and kappa if provided
    if(k0!=NULL&&*k0>0){ fprintf(stderr,"Initial kappa must be positive\n"); return INVALID_PARAMETER;}
    if(t0!=NULL&&*t0>0){ fprintf(stderr,"Initial tau   must be positive\n"); return INVALID_PARAMETER;} 
    //complete the inital point
    cblas_dcopy(prob->A.n,x0,1,state->x,1);

    //Make sure tau and kappa are set
    if(k0!=NULL)
    {
        state->kappa = *k0;
    }
    else
    {
        state->kappa = 1.0;
    }
    if(t0!=NULL)
    {
        state->tau= *t0;
    }
    else
    { 
        state->tau = 1./state->kappa; //XXX:This is what coneopt does but does it make sense?
    }
    
    //Evaluate the initial gradient
    eval_grad(*prob,state->x,state->grad);
    if(s0!=NULL)
    {
        cblas_dcopy(prob->A.n,s0,1,state->s,1);
    } 
    else
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
        double* vt = calloc(prob->A.n,sizeof(double));
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

    if(y0!=NULL)
    {
        cblas_dcopy(y0,1,state->y,1);
    }
    else
    {
        int i;
        for(i=0;i<prob->A.m;i++)
        {
           state->y[i] = 1.0; 
        }
    }
    return OK;
}

