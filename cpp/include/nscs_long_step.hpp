#ifndef H_NSCS
#define H_NSCS
    #ifdef __cplusplus
    extern "C" {
    #endif

    #define csi long int // Toggle 64 bit int

//Structure to hold the problem definition
typedef struct
{
   //Problem size and cone definitions
   csi m;
   csi n; //Populated during initialization
   //Cone sizes 
   csi n_free; //Free variables
   csi n_pos;  //Positive variables
   csi n_soc_cones; //Number of second order cones
   csi* soc_cones; //Array of dimensions of second order cones
   csi n_sdp_cones; //Number of sdp variables
   csi* sdp_cones; //Array of sdp variable sizes
   csi n_exp_cones; //Number of exponential cones

   double Av;
   csi*   Ai;
   csi*   Ap;
   csi    nnzA;
    
   double*   b;
   double*   c;

   double nu;     //Complexity: populated during initialization

} problem_t; 

//Parameters structure
typedef struct
{
    int max_iter;                //Maximum number of iterations
    int max_affine_backtrack_iter; //Maximum number of backtracking steps
    double backtrack_affine_constant;   //Affine backtracking constant
    double eta;         //Backtrack from the boundary

    double stop_primal;  //Stopping criteria bounds
    double stop_dual;
    double stop_gap;
    double stop_mu;
    double stop_tau_kappa;

    //True: uses the Mehrotra second order term
    bool solve_second_order;  
    //Degree of verbosity
    int  print;          
    // Regularization for the linear solver
    double delta;
    double gamma;
    int max_iter_ref_rounds; //Maximum number of iterative refinement 

} params_t;

typedef struct
{   
    double* xc;
    double* xf;
    double* y;
    double* s;

} state_t;
//ERRORS 
#define INVALID_PAR -8

    #ifdef __cplusplus
                }
    #endif
//Non C callable functions
void allocate_state(state_t &st,problem_t pr);
void print(params_t pars, int level, string s);
pars_t set_default_pars_nscs_long_step();

#endif
