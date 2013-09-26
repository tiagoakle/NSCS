#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "check.h"
#include "test_centering_measure.h"
#include "eval_cent_meas.h"
#include "OpenBLAS/cblas.h"
#include "test_util.h"
#include "smatvec.h"
#include "umfpack.h"
#include "spmat.h"
#include "linear_solvers.h"
#include "barriers.h"

START_TEST (test_centering_measure)
{

    //The centering measure squared is defined as 
    //||s-mug(x)||^2_H(x)^{-1}
    //Which is equivalent to s'H^{-1}s -2mu x's + mu^2 nu
    //if we choose s = H(x)\hat s 
    //This becomes s'\hat s -2 mu x's + mu^2 ny 
    
    //We test the centering mesure by generating a sequence of 
    //random feasible points for a problem with positive constraints
    //and exponential cones.
    //We calculate the centrality in the alternative way and compare
    
    //Define the problem structure
    problem_t prob;
    spmat A; 
    
    int n  = 16;
    A.n    = n;
    prob.A = A;
    int tk[] = {0,3,3};
    prob.tK = tk;
    csi nK[] = {10,3,3};
    prob.nK = nK;
    prob.k_count = 3;
    
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
    
    double *x = (double*)calloc(n,sizeof(double));
    double *s = (double*)calloc(n,sizeof(double));
    double *shat = (double*)calloc(n,sizeof(double)); 
    
    //Create the working vectors for the centrality measure
    double *psi =  (double*)calloc(n,sizeof(double)); 
    double *hpsi = (double*)calloc(n,sizeof(double)); 

    
    //Allocate the space to compress to csr
    int* Hi=(int*)calloc(nnzH,sizeof(int));
    int* Hp=(int*)calloc(n+1,sizeof(int));
    double* Hx=(double*)calloc(nnzH,sizeof(double)); 
    int* Map=NULL;
    
    //Define other vars
    int status;
    int i,j;
    double sshat;
    double mu;
    double nu = 16;

    double cent, cent_alt, cent_dif;
    //How many tests 
    int tests = 100;

    for(j=0;j<tests;j++)
    {
        //Pick a random s, 
        for(i=0;i<10;i++) s[i] = rand()/(double)RAND_MAX; 
        
        s[i] = rand()/(double)RAND_MAX - 0.5;
        s[i+2] = rand()/(double)RAND_MAX+1.;
        s[i+1] = -s[i+2]*exp(s[i]/s[i+2]-1) + rand()/(double)RAND_MAX;
        i = i+3;
 
        s[i] = rand()/(double)RAND_MAX - 0.5;
        s[i+2] = rand()/(double)RAND_MAX+1.;
        s[i+1] = -s[i+2]*exp(s[i]/s[i+2]-1) + rand()/(double)RAND_MAX;
       
        
        //Generate a random feasible point 
        for(i=0;i<10;i++) x[i] = rand()/(double)RAND_MAX; 
        
        x[i] = rand()/(double)RAND_MAX - 0.5;
        x[i+2] = rand()/(double)RAND_MAX+1.;
        x[i+1] = x[i+2]*exp(x[i]/x[i+2]) + rand()/(double)RAND_MAX;
        i = i+3;
 
        x[i] = rand()/(double)RAND_MAX - 0.5;
        x[i+2] = rand()/(double)RAND_MAX+1.;
        x[i+1] = x[i+2]*exp(x[i]/x[i+2]) + rand()/(double)RAND_MAX;
        
        //evaluate the hessian
        eval_hess(prob,x,state);
      
        //Calculate s 
        solve_linear_system(shat, H.I, H.J, H.V, nnzH, s , n);
        
        sshat = cblas_ddot(n,s,1,shat,1);
        mu    = cblas_ddot(n,s,1,x,1);
        mu = mu/nu;
        cent = sqrt(sshat -  mu*mu*nu); 
        cent_alt = eval_cent_meas(prob,x,s,state,mu,psi,hpsi);
        cent_dif  = fmax(cent-cent_alt,cent_alt-cent);

        ck_assert_msg(cent_dif<1.e-12,"In test %i centrality measures do not match |c-c_alt| %g\n",j,cent_dif);

    }

  
}
END_TEST

Suite* centering_measure_suite(void)
{
    Suite* suite = suite_create("Centering Measure");
    TCase *tc = tcase_create("Case1");
    tcase_add_test(tc,test_centering_measure);
    suite_add_tcase(suite,tc);
    return suite;
}


