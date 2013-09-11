#include "nscs.h"
#include "nscs_sp.h"
#include <math.h>
#include "barriers.h"
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
 * Where 
 * each cone is one of: positive orthant (0), second order cone (1), 
 * semi-definite cone (2), exponential cone (3) , power cone (4).
 * 
 * tK is an array of integers and of length q. The value tK[0] indicates the type
 * of cone 1,..., tK[q-1] indicates the type of cone K_q.
 *
 * nK is a length q array of integers. It parametrizes the size of the cones,
 * nK[0] the size of cone 1,..., nK[q] the size 
 * of cone q.
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
    
    //Validate the parameters
    int status = validate_pars(problem,pars);
    if(status != OK)
    {
        return status;
    }

    //This variable holds a reference to the state structure
    state_t state;
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
        eval_hess(*problem,state.x->v,state);
        //Calculate the approximate direction
        res = solve_approximate_tangent_direction(state);
        //Do a line search 
        res = linesearch_atd(state,*pars,*problem);
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
   
    if(!((state.H.I==NULL)&&(state.H.J==NULL)&&(state.H.V==NULL))) free_spmat(state.H);

}

/**
 * Validates the parameters provided by the user
 */
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


//TODO: THIS IS INCOMPLETE
void init(problem_t prob, state_t state, parameters_t params)
{
     //Count the number of non zeros in the new hessian
    int k = 0;
    csi nnzH = 0; 
    for(k=0;k<prob.k_count;k++)
    {
        nnzH = cone_nnz(prob.tK[k], prob.nK[k]);
    }
    //Allocate the new hessian
    //state.H = calloc_spmat(problem->n,problem->n,nnzH);
}

