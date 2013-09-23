/**
 * This file defines the structures
 * for nscs except for those which
 * define the sparse matrices and the vectors.
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

   //Present residuals
   double*  p_res; //Primal residual Ax-tb=p_res TODO:Check sign
   double*  d_res; //Dual residual -A'y+s-tc TODO:Check sign
   double*  c_res; //Complementarity res TODO:Define
   double g_res; //b'y-c'x-k = g_res TODO:Check sign
   double c_cent_res; // t/k - mu = c_rent_res TODO: Check this definition

   //Present residual norms
   double n_p_res;
   double n_d_res;
   double n_c_res;

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

   //Present value of the complementarity
   double mu;
 
   //Minimum centrality value found
   double min_centmeas;

   //Counters
   int nbacktrack;

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
   double  delta; //Primal regularization
   double  gamma; //Dual regularization
   
   double nu;     //Complexity
   csi    nnzH;   //Number of non zeros in the hessian

} problem_t; 

//Enumeration of problem types for the result strucutre
enum status_e
{
    optimal,
    p_infeasible,
    d_infeasible,
    mal_formed
};

//Result structure definition
typedef struct
{
  double* x;
  double* y;
  double* s;
  double tau;
  double kappa;

  enum status_e status;
  char*    status_message;
} result_t;

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

} parameters_t;

int nscs(problem_t* problem, parameters_t* pars, result_t* result);
void free_state(state_t state);
int validate_pars(problem_t* problem, parameters_t* pars);
int allocate_state(state_t state,problem_t* problem);
int  calculate_initial_point(state_t state,parameters_t* pars);
int  solve_approximate_tangent_direction(state_t state);

bool check_centering_condition(state_t state, parameters_t* pars);
int calculate_centering_direction(state_t state, parameters_t* pars);
int centering_linesearch(state_t state, parameters_t* pars);
void print_and_log(state_t state ,int m_iter, int c_iter,parameters_t* pars);
int check_stopping_criteria(state_t state,parameters_t* pars);
void print_final(state_t state);
void build_result(state_t state ,result_t* res, parameters_t* pars);

int calculate_residuals(state_t state);
#endif
