//Tests for the methods defined in nscs.c
//Test the functions that evaluate the barriers
#include "barriers.h"
#include "check.h"
#include "stdio.h"
#include "OpenBlas/cblas.h"
#include "test_nscs.h"
#include "nscs.h"
#include "smatvec.h"
#undef I


START_TEST(test_validate_range)
{
    double p,l,u;
    p = 1.;
    l = 0.5;
    u = 1.5;
    char* msg = "p";
    //In range
    ck_assert_msg(val_range(p,l,u,msg)==OK,"In range failed");
   //At lower
    p = l;
    ck_assert_msg(val_range(p,l,u,msg)==OK,"At lower failed");
   //At upper 
    p = u;
    ck_assert_msg(val_range(p,l,u,msg)==OK,"At upper failed");
   //Below lower
    p = l-1.e-10;
    ck_assert_msg(val_range(p,l,u,msg)!=OK,"Below lower failed ");
    //Above upper
    p = u+1.e-10;
    ck_assert_msg(val_range(p,l,u,msg)!=OK,"Above upper failed");

}
END_TEST

START_TEST(test_init)
{
    //Test the generation of s in the init routine
    problem_t prob;
    spmat A; 
    parameters_t pars;
    pars.delta = 1.;
    pars.gamma = 1.;
    
    int n  = 16;
    int m  = 10;

    A.n    = n;
    A.m    = m;
    prob.A = A;
    int tk[] = {0,3,3};
    prob.tK = tk;
    csi nK[] = {10,3,3};
    prob.nK = nK;
    prob.k_count = 3;
    prob.nu      = 10+6;
    
    int nnzH = 28;
    
    spmat H;
    H.n = A.n;
    H.m = A.n;
    H.I = (csi*)calloc(nnzH,sizeof(csi));
    H.J = (csi*)calloc(nnzH,sizeof(csi));
    H.V = (double*)calloc(nnzH,sizeof(double));
    H.nnz = nnzH;
    
    state_t state;
    state.H = H;
    allocate_state(&state,&prob);
    
    double *x = (double*)calloc(n,sizeof(double));        
    double *s = (double*)calloc(n,sizeof(double));        
    double *y = (double*)calloc(m,sizeof(double));        
    double tau;
    double kappa;
   

    //Generate a random feasible point 
    int i;
    for(i=0;i<10;i++) x[i] = rand()/(double)RAND_MAX; 
    
    x[i] = rand()/(double)RAND_MAX - 0.5;
    x[i+2] = rand()/(double)RAND_MAX+1.;
    x[i+1] = x[i+2]*exp(x[i]/x[i+2]) + rand()/(double)RAND_MAX;
    i = i+3;
 
    x[i] = rand()/(double)RAND_MAX - 0.5;
    x[i+2] = rand()/(double)RAND_MAX+1.;
    x[i+1] = x[i+2]*exp(x[i]/x[i+2]) + rand()/(double)RAND_MAX;

    int status; 
    status = init(&prob, &state, &pars, y, x, &tau, s, &kappa, false, false, false, false);
    ck_assert_msg(status==OK,"Init returned error");
    //At this point s must satisfy s+mug(x) = 0
    double * res = (double*)calloc(n,sizeof(double));

    double mu = cblas_ddot(n,state.x,1,state.s,1);
    mu = mu/(prob.nu);
    eval_grad(prob,state.x,res);
    cblas_dscal(n,mu,res,1);
    cblas_daxpy(n,1.,state.s,1,res,1);
    double n_res = cblas_dnrm2(n,res,1);
  
    ck_assert_msg(n_res < 1.e-15,"Norm of ||s+mug(x)||>1.e-15, nres: %g",n_res);
    ck_assert_msg(dual_feas(prob,x),"Generated initial dual is not feasible");
    //int init(problem_t* prob, state_t *state, double* y0, double* x0, double* t0, double* s0, double* k0)
}
END_TEST

START_TEST(test_calculate_residuals)
{
    //Make a dense matrix A of size 
    //m by n, a random b and c, x,y,s and
    //calculate the residuals using the function in
    //nscs and using blas
   
    int n  = 16;
    int m  = 10;

    problem_t prob;
    state_t state;
    
    allocate_state(&state,&prob);
    

    double *s = (double*)calloc(n,sizeof(double));        
    double *x = (double*)calloc(n,sizeof(double));        
    double *y = (double*)calloc(m,sizeof(double));
    double *b = (double*)calloc(m,sizeof(double));
    double *c = (double*)calloc(n,sizeof(double));
    
    //Generate a random primal point 
    int i,j;
    for(i=0;i<16;i++) x[i] = 1.; //rand()/(double)RAND_MAX; 
    //Generate a random dual point 
    for(i=0;i<16;i++) s[i] = rand()/(double)RAND_MAX; 
   
    //Generate a kappa and tau
    double kappa = rand()/(double)RAND_MAX;
    double tau   = rand()/(double)RAND_MAX;

    //Generate a random y
    for(i=0;i<m;i++){y[i] = rand()/(double)RAND_MAX - 0.5;}
    
    //Generate a random b
    for(i=0;i<m;i++){b[i] = rand()/(double)RAND_MAX - 0.5;}
   
    //Generate a random c
    for(i=0;i<n;i++){c[i] = 1.;} //= rand()/(double)RAND_MAX - 0.5;}

    //Allocate the sparse A
    spmat A;  
    A.n    = n;
    A.m    = m;
    A.nnz  = n*m;
    
    //Allocate a dense matrix A
    double* dA = (double*)calloc(m*n,sizeof(double));
    //Allocate space for a sparse matrix A
    csi* AI = (csi*)calloc(m*n,sizeof(csi));
    csi* AJ = (csi*)calloc(m*n,sizeof(csi));
    double* AV = (double*)calloc(m*n,sizeof(double));

    //Generate random entries for A
    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
        {
            dA[i+j*m] = rand()/(double)RAND_MAX - 0.5;
            AV[i+j*m] = dA[i+j*m];
            AI[i+j*m] = i;
            AJ[i+j*m] = j;
        }

    //Add b,c and A to prob
    prob.b = b;
    prob.c = c;
    A.I = AI;
    A.J = AJ;
    A.V = AV;
    prob.A = A;

    //Add x,s,y,tau,kappa to state
    state.x = x;
    state.s = s;
    state.y = y;
    state.tau = tau;
    state.kappa = kappa;

    //Calculate the residuals with blas
    double* dp_res = (double*)calloc(m,sizeof(double));
    double* dd_res = (double*)calloc(n,sizeof(double));
    //-b'y+c'x+k
    double  dg_res1 = -cblas_ddot(m,b,1,y,1) + cblas_ddot(n,c,1,x,1) + kappa;
     
    //-Ax+tb
    cblas_dcopy(m,prob.b,1,dp_res,1);
    cblas_dgemv(CblasColMajor,  CblasNoTrans,  m,  n,
		 -1., dA, m, x, 1,  tau,  dp_res, 1);

    //A^Ty+s-tc
    cblas_dcopy(n,prob.c,1,dd_res,1);
    cblas_dgemv(CblasColMajor,  CblasTrans,  m,  n,
		 1., dA, m, y, 1,  -tau,  dd_res, 1);
    cblas_daxpy(n,1.,s,1,dd_res,1);
    
    //Call the nscs calculate residuals method
    calculate_residuals(&state,&prob);
    //Compare the two calculations
    
    for(i=0;i<n;i++) ck_assert_msg(fmax(dd_res[i] - state.d_res[i],-dd_res[i] + state.d_res[i])<1.e-15,\
                                    "Dual residual not equal to blas based calculation: %g, i: %i",\
                                    fmax(dd_res[i] - state.d_res[i],-dd_res[i] + state.d_res[i]),i);
    for(i=0;i<m;i++) ck_assert_msg(fmax(dp_res[i] -state.p_res[i],-dp_res[i] + state.p_res[i])<1.e-15,\
                                    "Primal residual not equal to blas based calculation: %g, i: %i",\
                                    fmax(dp_res[i] -state.p_res[i],-dp_res[i] + state.p_res[i]),i);
    double a_err = fmax(state.g_res-dg_res1,-state.g_res+dg_res1);
    ck_assert_msg(a_err<1.e-15,"Gap residual not equal to blas calculation: abs(err): %g, %g, %g",a_err,state.g_res,dg_res1);

    //v.rD = (-pars.A'*v.y - v.s + pars.c*v.tau);
    //int init(problem_t* prob, state_t *state, double* y0, double* x0, double* t0, double* s0, double* k0)
    free(dd_res);
    free(dp_res);
    free(AV);
    free(AJ);
    free(AI);
    free(dA);
    free(c);
    free(b);
    free(y);
    free(x);
    free(s);

}
END_TEST

Suite* nscs_test_suite(void)
{

    Suite* suite = suite_create("nscs");
    TCase *tc = tcase_create("validate parameters");
    tcase_add_test(tc,test_validate_range);
    tcase_add_test(tc,test_calculate_residuals);
    tcase_add_test(tc,test_init);
    suite_add_tcase(suite,tc);
    return suite;
}
