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
 * nK is an array of integers and of size q which parametrizes the size of the cones,
 * nK[0] is the number of free variables, nK[1] the size of cone 1,..., nK[q-1] the size 
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
 * Backtracking linesearch for the aproximate tangent direction
 */
int linesearch_atd(state_t state ,parameters_t  pars, problem_t prob)
{
    double a0 = 1.; //Initial step length
    double a  = a0;
//  % set intial step length 
//    a0 = 1.0;
//    a  = a0;
//
//    if v.dkappa < 0
//        kapmax = -v.kappa/v.dkappa;
//        a = min(a,kapmax);
//    end
//    if v.dtau < 0
//        taumax = -v.tau/v.dtau;
//        a = min(a,taumax);
//    end
//    % if either kap or tau is blocking, multiply by eta
//    % so we do not hit boundary:
//    if a < a0
//        a = pars.eta*a;
//    end
    if(state.dkappa < 0) a = fmin(a,-state.kappa/state.dkappa);
    if(state.dtau < 0) a = fmin(a,-state.tau/state.dtau);
    if(a<a0) a = a*pars.eta;

//    % couter for number of bisections
//    nsect = 0;
//    
//    %Main linesearch loop
//    for j = 1:pars.lsmaxit


//  allocate work vectors for the trial steps
    double* xa = calloc(state.x->n,sizeof(double));
    double* sa = calloc(state.x->n,sizeof(double));
    double* psi  = calloc(state.x->n,sizeof(double));
    double* hpsi = calloc(state.x->n,sizeof(double));
    double kappaa;
    double taua;
    double dga;
    double mua;

    //TODO: add clean up before return?
    if(xa == NULL)   return OUT_OF_MEMORY;
    if(sa == NULL)   return OUT_OF_MEMORY;
    if(psi == NULL)  return OUT_OF_MEMORY;
    if(hpsi == NULL) return OUT_OF_MEMORY;

    int nsect = 0;
    int j = 0;

    //Define some variables
    //used in the loop.
    bool dFeas, pFeas; //Flags to indicate feasiblity
    bool dosect;       //Flag to indicate if there should be a backtrack

    double centmeas;   //Present value of centrality measure

    for(j=0;j<pars.max_backtrack;j++)
    {
//        xa     = v.x     + a * v.dx;
//        sa     = v.s     + a * v.ds;
//        taua   = v.tau   + a * v.dtau;
//        kappaa = v.kappa + a * v.dkappa;
      cblas_dcopy(state.x->n,state.x->v,1,xa,1);
      cblas_dcopy(state.s->n,state.s->v,1,sa,1);

      cblas_daxpy(prob.n,a,state.dx->v,1,xa,1);
      cblas_daxpy(prob.m,a,state.ds->v,1,sa,1);
      taua      = state.tau + a*state.dtau;
      kappaa    = state.kappa + a*state.dkappa;
        
//        % new duality gap:
//        dga    = xa'*sa + taua*kappaa;
//        mua    = dga / (K.nu + 1);
      dga    = cblas_ddot(state.x->n,xa,1,sa,1) + taua*kappaa;
      mua    = dga / (prob.nu + 1);  //TODO:Populate prob.nu on init

    //TODO:
    //Move the calculation of the nnzH to the initialization 
    //so that problem.nnzH is populated correctly
    //Add the calculation of problem.nu 
    //Add the population of problem m,n;
    //
    //Check if x,s are feasible wrt the cones
    dFeas = dual_feas(prob,state.s->v);
    pFeas = primal_feas(prob,state.x->v);

//        % evaluate barriers at new point:
//        % check only feasibility, so want = [-1,-1,-1]:
//        % Evaluate the primal barrier for f,g,H
//        FP = BarrFuncP(xa,K,[1,1,1]); 
//        % Check the dual feasibility
//        FD = BarrFuncD(sa,K,[1,-1,-1]);
//        
//        dosect = false; %True if we must backtrack
//        %If either the primal is infeasible or 
//        % the dual is infeasible backtrack
//        if FP{4} < 0 
//            dosect  = true;
//            R.block = 'pf';
//        elseif FD{4} < 0 
//            dosect  = true;
//            R.block = 'df';
//        else %If the iterate is pirmal and dual feasible evaluate the centrality
//            psi       = sa + mua*FP{2};
//            centmeas5 = sqrt(psi'*(FP{3}\psi)); %XXX: Linear solve
//            centmeas = centmeas5;
//            
//            if centmeas > mua*pars.theta
//                dosect  = true;
//                R.block = 'ce';
//            end
//        end

    if(!dFeas){ dosect = true;}
    else if(!pFeas){dosect = true;}
    else
    {
        //Evaluate the gradient at the preset point
        eval_grad(prob,xa,psi);
        //scale by mu and add s to psi
        cblas_dscal(prob.n,mua,xa,1);
        cblas_daxpy(prob.n,1.0,sa,1,psi,1);
        
        //Evaluate the hessian at the present test point
        eval_hess(prob,xa,state);
        
        //Solve the linear system H(hpsi)=psi
        int ret = solve_linear_system(hpsi,state.H.I,state.H.J,state.H.V,state.H.nnz,psi,state.H.n);
        if(ret!=0) return INTERNAL_ERROR; //XXX:We should check for numerical error and not crash 
        centmeas = cblas_ddot(psi,1,hpsi,1);
        
        //Decide if we need to backtrack
        if(centmeas > mua*pars.theta)
        {
            dosect = true;
        }
        
//        
//        if dosect
//            a     = a*pars.lscaff; 
//            nsect = nsect + 1;
//        else
//            break;
//        end
//        
//    end %end of main linesearch loop
//

    }

    if(dosect)
    {
        a = a*pars.lscaff;
        state.nbacktrack += 1; 
    }
     

}

    //Free the work vectors 
    free(xa);
    free(sa);
    free(psi);
    free(hpsi);

//    v.a = a;
//    
//    % Check if the linesearch did not find 
//    % a feasible point in the maximum number of iterations
//    if j == pars.lsmaxit
//        xa = v.x + a * v.dx;
//        FP = BarrFuncP(xa,K,[1,1,-1]);
//        if FP{4} < 0
//
//            error(['linesearch: failed to find feasible point.',...
//                ' pars.lscaff too close to 1 ???']);
//        end
//    end
//
//    %Take the step 
//
//    % store the previous step
//    v.xprev     = v.x;
//    v.gprev     = v.F{2};
//    v.tauprev   = v.tau;
//    v.kappaprev = v.kappa;
//    
//    % take step:
//    v.x     = xa;
//    v.tau   = taua;
//    v.y     = v.y     + v.a * v.dy;
//    v.s     = sa;
//    v.kappa = kappaa;
//    
//    % update other quantities:
//    v.dgap  = v.x'*v.s + v.tau*v.kappa;
//    v.mu    = v.dgap / (K.nu + 1);
//    
//    % "feas" measure, see Sturm:
//    v.feas  = v.dtauaff/v.tauprev - v.dkappaaff/v.kappaprev;
//
//    %Update the residuals
//    bty  = pars.b'*v.y;
//    ctx  = pars.c'*v.x;    
//    v.rA = abs( ctx - bty )/( v.tau + abs(bty) );
//    
//    v.rP = (pars.A*v.x - pars.b*v.tau);
//    v.rD = (-pars.A'*v.y - v.s + pars.c*v.tau);
//    v.rG = (-ctx + bty - v.kappa);
//    
//    v.rPrel = norm( v.rP, 'inf')/pars.relstopP;
//    v.rDrel = norm( v.rD, 'inf')/pars.relstopD;
//    v.rGrel = norm( v.rG, 'inf')/pars.relstopG;
//    v.rArel = v.rA; 
   
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

