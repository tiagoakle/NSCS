//Tests for the methods defined in nscs.c
//Test the functions that evaluate the barriers
#include "barriers.h"
#include "check.h"
#include "stdio.h"
#include "OpenBlas/cblas.h"
#include "test_nscs.h"
#include "nscs.h"
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
    H.I = calloc(nnzH,sizeof(csi));
    H.J = calloc(nnzH,sizeof(csi));
    H.V = calloc(nnzH,sizeof(double));
    H.nnz = nnzH;
    
    state_t state;
    state.H = H;
    allocate_state(&state,&prob);
    
    double *x = calloc(n,sizeof(double));        
   

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
    status = init(&prob, &state, NULL, x, NULL, NULL, NULL);
    ck_assert_msg(status==OK,"Init returned error");
    //At this point s must satisfy s+mug(x) = 0
    double * res = calloc(n,sizeof(double));

    double mu = cblas_ddot(n,state.x,1,state.s,1);
    mu = mu/prob.nu;
    eval_grad(prob,state.x,res);
    cblas_dscal(n,mu,res,1);
    cblas_daxpy(n,1.,state.s,1,res,1);
    double n_res = sqrt(cblas_ddot(n,res,1,res,1));
  
    ck_assert_msg(n_res < 1.e-15,"Norm of ||s+mug(x)||>1.e-15, nres: %g",n_res);
    ck_assert_msg(dual_feas(prob,x),"Generated initial dual is not feasible");
    //int init(problem_t* prob, state_t *state, double* y0, double* x0, double* t0, double* s0, double* k0)
}
END_TEST

Suite* nscs_test_suite(void)
{

    Suite* suite = suite_create("nscs");
    TCase *tc = tcase_create("validate parameters");
    tcase_add_test(tc,test_validate_range);
    tcase_add_test(tc,test_init);
    suite_add_tcase(suite,tc);
    return suite;
}
