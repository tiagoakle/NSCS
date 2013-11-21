/* Main file for the solver nscs_long_step
 */
#include "nscs_long_step.hpp"
#include <String>
using namespace std;

int nscs_long_step(problem_t p, double* x0c, double* x0f, params_t par)
{
    if(x0c==0)
    {
        print(par,1,"call to nscs_long_step with null initial x0c point");
        return INVALID_PAR;
    }
    //X0f can be null if all variables are constrained
    if(x0f==0&&p.n_free!=0)
    {
        print(par,1,"call to nscs_long_step problem_t.n_pos !=0 and x0f null");
        return INVALID_PAR;
    }
    if(par==0)
    {
        //Set the default parameters
        par = set_default_pars_nscs_long_step(); 
    }

    //Validate the problem structure
    int ret = validate_problem_structure(problem_t pr);
    if(ret != RET_OK)
    {
        print(par,"call to nscs_long_step with invalid problem structure");
        return INVALID_PAR;
    }

    //Initialize the state structure
    state_t st;
    allocate_state(&st,p);  
    return 0;
}

void allocate_state(state_t &st,problem_t pr)
{
}

void print(params_t pars, int level, string s)
{ 
    if(pars.print>level)
        printf(s);
}

pars_t set_default_pars_nscs_long_step()
{
    pars_t pars;
    pars.max_iter   = 100;  //Maximum outer iterations
    pars.max_affine_backtrack_iter = 300;  //Maximum affine backtracking steps
    pars.backtrack_affine_constant = 0.9;  //Affine backtracking constant
    pars.eta        = 0.98;                //Multiple of step to the boundary
    
    pars.stop_primal= 1e-5;   //Stopping criteria bounds
    pars.stop_dual  = 1e-5;
    pars.stop_gap   = 1e-5;
    pars.stop_mu    = 1e-7;
    pars.stop_tau_kappa = 1.e-7;

    pars.solve_second_order = true;

    pars.print      = 1;             //Level of verbosity from 0 to 11
    //Regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 20;

}
