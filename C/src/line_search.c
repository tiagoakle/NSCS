#include "nscs.h"
#include "nscs_sp.h"
#include <math.h>
#include "barriers.h"
#include "line_search.h"
#include "stdio.h"
#include "OpenBlAS/cblas.h"
#include "linear_solvers.h"

/**
 * Wrapper to call from matlab
 */
int linesearch_atd_no_structs( int m, int n, double*x, double*y, double*s, double tau, double kappa, double* dx,\
                                                                                                    double* dy,\
                                                                                                    double* ds,\
                                                                                                    double  dtau,\
                                                                                                    double  dkappa,\
                                                                                                    double  lscaff,\
                                                                                                    double  eta,\
                                                                                                    double  theta,\
                                                                                                    int max_backtrack,\
                                                                                                    int k_count,\
                                                                                                    int* nK,\
                                                                                                    int* tK,\
                                                                                                    double nu,\
                                                                                                    csi nnzH,\
                                                                                                    double* a,\
                                                                                                    int* nbacktrack)
{

    //Construct prob
    problem_t prob;
    prob.tK      = tK;
    prob.nK      = nK;
    prob.k_count = k_count;
    prob.m       = m;
    prob.n       = n;
    prob.delta   = 1.e-10;
    prob.gamma   = 1.e-10;
    prob.free    = 0;
    prob.nu      = nu;
    prob.nnzH    = nnzH;

    //Construct the state
    state_t state;
    vec vy;
    state.y = &vy;
    state.y->v = y;
    state.y->n = m;

    vec vx;
    state.x = &vx;
    state.x->v = x;
    state.x->n = n;
    
    vec vs;
    state.s = &vs;
    state.s->v = s;
    state.s->n = n;
    
    state.tau = tau;
    state.kappa = kappa;


    vec vdy;
    state.dy = &vdy;
    state.dy->v = dy;
    state.dy->n = m;

    vec vdx;
    state.dx = &vdx;
    state.dx->v = dx;
    state.dx->n = n;
    
    vec vds;
    state.ds = &vds;
    state.ds->v = ds;
    state.ds->n = n;
    
    state.dtau = dtau;
    state.dkappa = dkappa;
    
    state.mu     = NAN; //In case we access this
    state.eta    = NAN;
    state.nbacktrack = 0;
    
    //Build the params structure

    parameters_t params;
    params.print = true;
    params.max_center_iter = NAN;
    params.eta   = eta;
    params.theta = theta;
    params.lscaff = lscaff;
    params.max_backtrack = max_backtrack;

    //Allocate space for the hessian!!!
    state.H.I = (int*)calloc(nnzH,sizeof(int));
    state.H.J = (int*)calloc(nnzH,sizeof(int));
    state.H.V = (double*)calloc(nnzH,sizeof(double));

    //Call the linesearch
    linesearch_atd(state,params,prob);
    (*nbacktrack) = state.nbacktrack;
    (*a)          = state.a;
    
    free(state.H.I);
    free(state.H.J);
    free(state.H.V);

}

void print_max_abs(double*x, int n, char* name)
{
    double max = 0;
    int i;
    for(i=0;i<n;i++)
    {
        max = fmax(fmax(x[i],-x[i]),max);
    }
    printf("%s: %g\n",name,max);
}

void print_min(double*x, int n, char* name)
{
    double min = 0;
    int i;
    for(i=0;i<n;i++)
    {
        min = fmin(x[i],min);
    }
    printf("%s: %g\n",name,min);

}
/**
 * Backtracking linesearch for the aproximate tangent direction
 */
int linesearch_atd(state_t state ,parameters_t  pars, problem_t prob)
{
    printf("linesearch \n");
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
    double* xa = (double*)calloc(state.x->n,sizeof(double));
    double* sa = (double*)calloc(state.x->n,sizeof(double));
    double* psi  = (double*)calloc(state.x->n,sizeof(double));
    double* hpsi = (double*)calloc(state.x->n,sizeof(double));
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
         printf("Iterate %i,%g,",j,a);
//        xa     = v.x     + a * v.dx;
//        sa     = v.s     + a * v.ds;
//        taua   = v.tau   + a * v.dtau;
//        kappaa = v.kappa + a * v.dkappa;
      cblas_dcopy(state.x->n,state.x->v,1,xa,1);
      cblas_dcopy(state.s->n,state.s->v,1,sa,1);

      cblas_daxpy(prob.n,a,state.dx->v,1,xa,1);
      cblas_daxpy(prob.n,a,state.ds->v,1,sa,1);
      taua      = state.tau + a*state.dtau;
      kappaa    = state.kappa + a*state.dkappa;
        
//        % new duality gap:
//        dga    = xa'*sa + taua*kappaa;
//        mua    = dga / (K.nu + 1);
      dga    = cblas_ddot(state.x->n,xa,1,sa,1) + taua*kappaa;
      mua    = dga / (prob.nu + 1);  //TODO:Populate prob.nu on init
      
   //   //XXX:REMOVE
   //   if(j==0)
   //   {
   //   write_vector_to_csv("first_iter_sa.csv",sa,state.s->n);
   //   write_vector_to_csv("first_iter_ds.csv",state.ds->v,state.s->n);
   //   write_vector_to_csv("first_iter_s.csv" ,state.s->v,state.s->n);
   //   }

    //TODO:
    //Move the calculation of the nnzH to the initialization 
    //so that problem.nnzH is populated correctly
    //Add the calculation of problem.nu 
    //Add the population of problem m,n;
    //
    //Check if x,s are feasible wrt the cones
    dFeas = dual_feas(prob,sa);
    pFeas = primal_feas(prob,xa);

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

    
    if(!pFeas)
    {
        dosect = true; printf("C: Not primal feasible\n");
    }
    else if(!dFeas)
    {
        dosect = true; printf("C: Not dual feasible\n");
    }
    else
    {

        printf("Evaluating dist\n");
        //Evaluate the gradient at the preset point
        eval_grad(prob,xa,psi);
        //scale by mu and add s to psi
        cblas_dscal(prob.n,mua,xa,1);
        cblas_daxpy(prob.n,1.0,sa,1,psi,1);
         
        //Evaluate the hessian at the present test point
        eval_hess(prob,xa,state);
         
        fflush(stdout);
        //Solve the linear system H(hpsi)=psi
        int ret = solve_linear_system(hpsi,state.H.I,state.H.J,state.H.V,state.H.nnz,psi,state.H.n);
        if(ret!=0) return INTERNAL_ERROR; //XXX:We should check for numerical error and not crash 
        centmeas = sqrt(cblas_ddot(prob.n,psi,1,hpsi,1)); 

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
    state.a = a;   
    return OK;
}


