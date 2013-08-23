/**
 * This file defines the structures
 * for nscs except for those which
 * define the sparse matrices and the vectors.
 *
 */

#ifndef H_NSCS
#define H_NSCS
#include "common.h"
#include "nscs_sp.h"
//State structure
typedef struct
{

   //Present iterate 
   vec*   y;
   vec*   x;
   double tau;
   vec*   s;
   double kappa;

   //Present residuals
   vec*  p_res; //Primal residual Ax-tb=p_res TODO:Check sign
   vec*  d_res; //Dual residual -A'y+s-tc TODO:Check sign
   vec*  c_res; //Complementarity res TODO:Define
   double g_res; //b'y-c'x-k = g_res TODO:Check sign
   double c_cent_res; // t/k - mu = c_rent_res TODO: Check this definition

   //Present residual norms
   double n_p_res;
   double n_d_res;
   double n_c_res;

   //Present search direction
   vec*  dy;
   vec*  dx;
   double dt;
   vec*  ds;
   double dk;

   //Present scaling point
   vec* scal; 

   //Present value of the Hessian
   spmat* H; 

   //Present value of the complementarity
   double mu;
   //Present value of the centrality
   double eta;

   //This vector holds pointers to
   //the functions that return the number of 
   //non zeros for a Hessians of each cone type.
   //hessians_nnz[i](dim) returns the number of 
   //non zeros for a cone of type i of dimension dim 
   int (*hessians_nnz[CONE_TYPES])(int);
   
   //This vector holds pointers to 
   //the functions that calculate the hessian 
   //evaluated at the present scaling point.
   // The call 
   // nnz = (hessians[i])(first_ix,dim,state); 
   // returns the number of non zeros and
   // assigns the triplets I,J,V to state->H
   // starting at index first_ix. 
   int (*hessians[CONE_TYPES])(int,int,vec*);

} state_t;

//Problem definition
typedef struct
{
   //Problem definition
   spmat A;
   vec   b;
   vec   c;
   csi   *iK;
   int   *nK;
   int   m;
   int   n;
   int   free; //This is equal to iK[0] if nK[0] = 1
               // if nK[0] = 0 free = 0;
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
  vec x;
  vec y;
  vec s;
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
} parameters_t;

int nscs(problem_t* problem, parameters_t* pars, result_t* result);
void free_state(state_t state);
int validate_pars(problem_t* problem, parameters_t* pars);
int allocate_state(state_t state,problem_t* problem);
int  calculate_initial_point(state_t state,parameters_t* pars);
int eval_barrier(state_t state, problem_t* problem);
int  solve_approximate_tangent_direction(state_t state);


bool check_centering_condition(state_t state, parameters_t* pars);
int calculate_centering_direction(state_t state, parameters_t* pars);
int centering_linesearch(state_t state, parameters_t* pars);
int primal_lifting(state_t state);
void print_and_log(state_t state ,int m_iter, int c_iter,parameters_t* pars);
int check_stopping_criteria(state_t state,parameters_t* pars);
void print_final(state_t state);
void build_result(state_t state ,result_t* res, parameters_t* pars);

int linesearch(state_t state ,parameters_t * pars);
#endif
