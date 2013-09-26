/**
 * This file defines the structures
 * for nscs except for those which
 * define the sparse matrices 
 *
 */
#ifndef H_NSCS
#define H_NSCS
#include "common.h"
#include "spmat.h"
#include "math.h"
//State structure
typedef struct
{

   //Present iterate 
   double*   y;
   double*   x;
   double tau;
   double*   s;
   double kappa;

   //Previous iterates
   double prev_tau;
   double prev_kappa;

   //Present residuals
   double*  p_res; //Primal residual -Ax+tb=p_res 
   double*  d_res; //Dual residual +A'y+s-tc 
   double   g_res; //-b'y+c'x+k = g_res

   //Relative gap
   double rel_gap; //|c'x-b'y|/(tau+|b'y|)

   //Relative stop parameters
   double p_relstop;
   double d_relstop;
   double g_relstop;

   //Present relative residual norms they are of the form
   //||p_res||_inf/relstopP
   double n_p_res;
   double n_d_res;
   double n_g_res;

   //Present search direction
   double*  dy;
   double*  dx;
   double dtau;
   double*  ds;
   double dkappa;

   //Present step length
   double a;

   //Present value of the Hessian
   spmat H; 

   //Present value of the primal gradient
   double* grad;

   //Present value of the complementarity
   double mu;
   double dgap;
   //Feasibility measure! TODO:Justify
   double vfeas;

   //Centering dxnorm ||dx||_x
   double dx_norm;
 
   //Minimum centrality value found
   double min_centmeas;

   //Counters
   int m_iter;
   int nbacktrack;
   int ncentbacktrack;
   int ncent_iter;

   //Initial centrality
   double mu0;

} state_t;

//Structure to hold the problem definition
typedef struct
{
   //Constraint matrix
   spmat A;
   double*   b;
   double*   c;

   int   *tK;
   csi   *nK;
   int   k_count; //Number of cones
   
   double nu;     //Complexity
   //These are copied from pars 
   double  delta; //Primal regularization
   double  gamma; //Dual regularization   

} problem_t; 

//Enumeration of problem types for the result strucutre
enum status_e
{
    optimal,
    p_infeasible,
    d_infeasible,
    mal_formed
};

//Define the pars structure
typedef struct 
{
    bool print;
    int max_iter;
    int max_center_iter;
    double theta; //Backtracking limit
    double lscaff; //Affine scaling backtracking linesearch constant
    double lsccent; //Centering backtracking linesearch constant
    double eta;    //Closeness to boundary of first iterate
    int    max_backtrack;
    //Regularization parameters
    double delta;
    double gamma;
    //Centering measure
    double beta;
    
    //Stopping parameters
    double p_relstop;
    double d_relstop;
    double rel_gap_relstop;
    double rho_g;
    double rhoI;
    double rhoM;

} parameters_t;

int nscs(problem_t* problem, parameters_t* pars,\
        double* y, double* x, double* t, double* s, double* k, bool wy, bool wt, bool ws, bool wk);
void free_state(state_t state);
int validate_pars(problem_t* problem, parameters_t* pars);
int allocate_state(state_t* state,problem_t* problem);
int  solve_approximate_tangent_direction(state_t *state,problem_t prob,parameters_t pars);
bool check_centering_condition(state_t *state, parameters_t* pars);
int calculate_centering_direction(state_t *state, problem_t prob, parameters_t* );
void print_and_log(state_t *state,parameters_t* pars);
void print_header();
void print_final(int stop_reason);
int check_stopping_criteria(state_t *state,parameters_t* pars,problem_t* prob);
void calculate_residuals(state_t* state, problem_t* prob);
int val_range(double param, double min, double max, char* name);
int init(problem_t* prob, state_t *state, parameters_t *pars, double* y0,\
        double* x0, double* t0, double* s0, double* k0, bool wy, bool wt, bool ws, bool wk);
#endif
