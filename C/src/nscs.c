#include "nscs.h"
#include "nscs_sp.h"

/**
 * Main entry point for the Non symmetric cone solver.
 *
 *
 * NSCS solves a conic programming problem given in 
 * standard form,
 * min c'x st Ax=b x_c \in K.
 * Where x_c denotes the set of constrained entries of x.
 * K is the product of several cones defined by the
 * arrays nK, iK. 
 * 
 * K = K_1xK_2xK_3,.....,K_q
 * 
 * Where each cone is one of: positive orthant, second order cone, 
 * semi-definite cone, exponential cone, power cone.
 *
 * nK indicates how many instances of each cone appear in the 
 * expression for K.
 *
 * nK[0] positive orthant cones can be 0 or 1.
 * nK[1] number of second order cones.
 * nK[2] number of semi-definite cones.
 * nK[3] number of exponential cones.
 * nK[4] number of power cones.
 *
 * The sum of nK must equal q
 *
 * iK is a (q+1) sized array which holds the cumulative count of variables
 * in the cones.
 * 0,...,nK[0]-1     indices of free varaibles
 * nK[0],...,nK[1]-1 indices of the variables in the positive orthant
 * nK[1],...,nK[2]-1 indices of variables in cone K_2
 * .
 * .
 * .
 * nK[q-1],...,nK[q]-1 indices of variables in cone K_q
 *
 * The variables with indices 0,..,K[0]-1 are free variables
 * K[0],...,K[1]-1 non negative
 * K[1],...,K[2]-1 belong to the second order cone
 * 
 *
 * @param A the constraint matrix
 * @param b the vector b 
 * @param c the vector c
 * @param K the cone definition array
 * @param pars msic prameters
 * @param res result
 *
 */
int nscs(problem_t* problem, parameters_t* pars, result_t* result) 
{

    //Initialize some stuff
    int res;

    //This variable holds a reference to the state
    state_t state;
    
    //Validate the parameters
    int status = validate_pars(problem,pars);
    if(status != OK)
    {
        return status;
    }

    //Allocate the structure to hold the state of the algorithm
    status    =  allocate_state(state, problem);
    if(status != OK)
    {
        free_state(state);
        return status;
    }

    //Calculate the starting iterate and save it in state
    res = calculate_initial_point(state, pars); 
    
    //Main loop --------------------------------------------
    int m_iter = 0; //Main iteration counter
    
    for(m_iter = 0;m_iter < pars->max_iter; m_iter++)
    {
        //Evaluate the barrier function
        res = eval_barrier(state,problem);
        //Calculate the approximate direction
        res = solve_approximate_tangent_direction(state);
        //Do a line search 
        res = linesearch(state,pars);
        //Calculate the new residuals
        res = calculate_residuals(state);
         
        bool do_center = !check_centering_condition(state,pars);
        int c_iter = 0;
        //Run the centering routine
        while(do_center && c_iter < pars->max_center_iter)
        {
            res = calculate_centering_direction(state,pars);
            //Linesearch for centering 
            res = centering_linesearch(state,pars);
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

        //At this point the iterate is centered and 
        //we need to do the lifting to calculate
        //the approximate scaling point
        res = primal_lifting(state);
        
        //Print the output and check for exit
        if(pars->print)
        { 
            print_and_log(state,m_iter,c_iter,pars);
        }
       
       //Check the stopping criteria
       res = check_stopping_criteria(state,pars);

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
    if(!(state.dy==NULL))   free_vec(state.dy);
    if(!(state.dx==NULL))   free_vec(state.dx);
    if(!(state.ds==NULL))   free_vec(state.ds);
   
    //(Free the residual vectors
    if(!(state.p_res==NULL))   free_vec(state.p_res);
    if(!(state.d_res==NULL))   free_vec(state.d_res);
    if(!(state.c_res==NULL))   free_vec(state.c_res);
   
    if(!(state.scal==NULL)) free_vec(state.scal);
    if(!(state.H==NULL))    free_spmat(state.H);
}


int validate_pars(problem_t* problem, parameters_t* pars)
{
    return VALIDATION_OK;
}

//Calculate the starting iterate and save it in state_t
int allocate_state(state_t state,problem_t* problem)
{
    //Allocate the state structure
    //state = calloc(1,sizeof(state_t));
    //if(state == NULL) return OUT_OF_MEMORY;
    
    //Allocate the vectors that hold the present iterate
    state.y = calloc_vec(problem->m);
    if(state.y == NULL) return OUT_OF_MEMORY;
    state.x = calloc_vec(problem->n);
    if(state.x == NULL) return OUT_OF_MEMORY;
    state.s = calloc_vec(problem->n-problem->free); //s is the size of the dimensions of all constrained variables
    if(state.s == NULL) return OUT_OF_MEMORY;

    //Allocate the vectors that hold the present residual
    state.p_res = calloc_vec(problem->m);
    if(state.p_res == NULL) return OUT_OF_MEMORY;
    state.d_res = calloc_vec(problem->n);
    if(state.d_res == NULL) return OUT_OF_MEMORY;
    state.c_res = calloc_vec(problem->n-problem->free);
    if(state.c_res == NULL) return OUT_OF_MEMORY;
    
    return OK;
}

int  calculate_initial_point(state_t state,parameters_t* pars)
{
    return 1;
}

int eval_barrier(state_t state, problem_t* problem)
{
    //Clear the present hessian
    if(state.H==NULL) free_spmat(state.H);
    
    //Build the new hessian
    int k = 0;
    int nnzH = 0; 
    //In this first pass we only calculate the nnz of the hessian
    //Iterate over all the cones and build the hessian for each 
    int k_count = 0; 
    int k_type  = 0; 
    for(k_type = 0;k_type<CONE_TYPES;k_type++)
    {
        //For each cone of each type
        for(k = 0;k<problem->nK[k_type];k++)
        {
            int (*foo)(int)  =  (state.hessians_nnz[k_type]); //TODO try to get rid of foo
            nnzH += foo(problem->iK[k_count]);  //TODO: Figure out the cast and call
            k_count++;
        }
    }
    
    //Allocate the space for the Hessian
    state.H = calloc_spmat(problem->n-problem->free,\
                           problem->n-problem->free,\
                           nnzH);

    //In this second pass generate the Hessian
    csi nnz    = 0; //How many non zeros have we generated
    
    //Iterate over all the cones and build the hessian for each 
    for(k_type = 0;k_type<CONE_TYPES;k_type++)
    {
        //For each cone of each type
        for(k = 0;k<problem->nK[k_type];k++)
        { 
            int (*foo)(int,int,vec*) = state.hessians[k_type]; //TODO: Try to get rid of foo
            nnz += foo(nnz,problem->iK[k_count],state.scal); //TODO: Figure out the cast and call
            k_count++;
        }
    }
 
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
 */
int linesearch(state_t state ,parameters_t * pars)
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

int primal_lifting(state_t state)
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
 
