//Test the functions that evaluate the barriers
#include "barriers.h"
#include "check.h"
#include "stdio.h"
#include "cblas.h"
#include "test_barriers.h"

START_TEST(s_test_pos_orthant)
{
    int n = 100;
    double* x = calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1+i*1.0;
    }
    ck_assert(pos_orthant_feas(x,n));
    
    x[1] = -1.0;
    ck_assert(!pos_orthant_feas(x,n));
    free(x);
}
END_TEST

//Test gradient for the positive orthant barrier 
START_TEST(s_test_pos_orthant_grad)
{
    int n = 100;
    double* x = calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1+i*1.0;
    }
    double* g = calloc(n,sizeof(double));
    pos_orthant_grad(g,x,n);
    
    for(i=0;i<n;i++)
    {
        if((g[i] - 1./x[i])>1.e-16) break;
    }
    
    printf("%i difference\n",i);
    ck_assert(i==n); 
    free(x);
    free(g);
}
END_TEST


//Test the hessian for the positive orthant barrier
START_TEST(s_test_pos_orthant_hess)
{
    int n = 100;
    double* x = calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1+i*1.0;
    }

    double delta = 1;

    csi* HI = calloc(n,sizeof(csi));
    csi* HJ = calloc(n,sizeof(csi));
    double* HV = calloc(n,sizeof(double));
    pos_orthant_hessian(HI,HJ,HV,x,n,delta);
    
    for(i=0;i<n;i++)
    {
        if((HV[i] - (1./x[i])*(1./x[i])-delta)>1.e-16) break;
    } 
    ck_assert_msg(i==n,"expected 1/x^2+d but entry %i is %g for x=%g\n",\
                        i,(1./x[i])*(1./x[i])+delta,x[i]);

    for(i=0;i<n;i++)
    {
        if(HI[i] != i) break;
    } 
    ck_assert_msg(i==n,"expected HI[i] = i but H[%i]=%i",\
                        i,HI[i]);
 
    for(i=0;i<n;i++)
    {
        if(HJ[i] != i) break;
    } 
    ck_assert_msg(i==n,"expected HJ[i] = i but H[%i]=%i",\
                        i,HJ[i]);


    free(x);
    free(HI);
    free(HJ);
    free(HV);
}
END_TEST


//Test the value of the barrier, evaluate at a known point
START_TEST(s_test_pos_orthant_val)
{

    int n = 100;
    double* x = calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1;
    }
    ck_assert(pos_orthant_val(x,n)==0);
    free(x);

}
END_TEST


//Test the nnz of the positive orthant hessian
START_TEST(s_test_pos_orthant_nnz)
{
    csi n = 100;
    ck_assert(pos_orthant_nnz(n)==n);
}
END_TEST

//Test the function that returns the complexity 
//parameter of the pos orthatnt barrier
START_TEST(s_test_pos_orthant_complexity)
{
    csi n = 100;
    ck_assert(pos_orthant_complexity(n)==n);
}
END_TEST

//Test the exponential barrier
START_TEST(test_exp_nnz)
{
    ck_assert(exp_nnz() == 9);
}
END_TEST

START_TEST(test_exp_complexity)
{
    ck_assert(exp_complexity() == 3);
}
END_TEST

START_TEST(test_exp_val)
{
    double x[] = {0.5,exp(1),1.};
    double val = exp_val(x);
    ck_assert(val-(-log(0.5)-1)<=1.e-15);
    ck_assert(-val+(-log(0.5)-1)<=1.e-15);
}
END_TEST

START_TEST(test_exp_grad)
{
    double x[] = {0.5,exp(1),1.};
    double xt[] = {0,0,0};
    double g[] = {0.,0.,0.};
    exp_grad(g,x);
    
    double s = 1.e-22;
    int j;
    double tm;
    
    for(j=0;j<3;j++) xt[j] = x[j]; //Copy x to the test vector        
    xt[0] = xt[0]+s; 
    tm = exp_val(xt)-exp_val(x);
    tm = tm/s-g[0];
    ck_assert_msg(tm<1.e-15,"Numerical approx g[0] %e\n",tm);
    
    double tm_min = 1;
    for(s = 1e-2;s>1.e-16;s=s/2)
    {
        for(j=0;j<3;j++) xt[j] = x[j]; //Copy x to the test vector        
        xt[1] = xt[1]+s; 
        tm = exp_val(xt)-exp_val(x);
        tm = tm/s-g[1];
        tm_min = fmin(tm,tm_min);
    }
    ck_assert_msg(tm_min<1.e-7,"Numerical approx g[1] %e\n",tm_min);

    tm_min = 1;
    for(s = 1e-2;s>1.e-16;s=s/2)
    {
        for(j=0;j<3;j++) xt[j] = x[j]; //Copy x to the test vector        
        xt[2] = xt[2]+s; 
        tm = exp_val(xt)-exp_val(x);
        tm = tm/s-g[2];
        tm_min = fmin(tm,tm_min);
    }
    ck_assert_msg(tm_min<1.e-7,"Numerical approx g[2] %e\n",tm_min);

}
END_TEST

START_TEST(test_exp_hess)
{
   csi HI[9];
   csi HJ[9];
   double HV[9];
   double x[] = {0.5,exp(1),1.};
   double xt[3];
   double g[3];
   double gt[3];
   double Hs[3];
   exp_hessian(HI,HJ,HV,x,0.); 

   //convert the hessian into compressed column form
 
     int Ai[9];
     int Ap[4];     
     double Av[9]; 
     int* Map=NULL;
     //Execute the call
     int status;
     status  =  umfpack_di_triplet_to_col(3,3,9,HI,HJ,HV,Ap,Ai,Av,Map);
     ck_assert_msg(status==0,"Unable to generate csr form");
   //Each test makes a product with the hessian and compares with the gradient
    int tests = 1;
    int t = 0;
    int i = 0;
    double a = 1; 
    double s[3];
    double n_sq;
    double s_norm;

    double max_error = 0.;
    double min_err_ap;
    for(t=0;t<tests;t++)
    {
        

        //Build a random direction 
        for(i=0;i<3;i++)
            s[i] = rand() / (double)RAND_MAX - 0.5;
        //Normalize
        s_norm = cblas_ddot(3,s,1,s,1);
        s_norm = sqrt(s_norm);
        cblas_dscal(3,1./s_norm,s,1);

        //Calculate the approximation
        min_err_ap = 1.; 
        for(a=1;a>1e-10;a=a/2){

            cblas_dcopy(3,x,1,xt,1); //copy x into xt
            cblas_daxpy(3,a,s,1,xt,1); //Shift xt by a*s

            //Evaluate the gradients
            exp_grad(g,x);
            exp_grad(gt,xt);
       
            //Test
            double sol[3];
            sol[0] = 0;
            sol[1] = 0;
            sol[2] = 0;
            Hs[0]  = 0;
            Hs[1]  = 0;
            Hs[2]  = 0;


            //Evaluate the product Hs
            dspmv(3,3,1.,Ap,Ai,Av,Hs,s); 
           
            //calculate -1/a(a*Hs+gs-gt)
            //as Hs+(-g+gt)(-1/a)
            cblas_daxpy(3,-1.,g,1,gt,1);
            cblas_daxpy(3,-1./a,gt,1,Hs,1);

            //Calculate the norm squared
            n_sq = cblas_ddot(3,Hs,1,Hs,1);
            min_err_ap = fmin(min_err_ap,sqrt(n_sq));//keep the smallst error

        }
        max_error = fmax(min_err_ap,max_error); //Find the largest over all the random tests
    }
    ck_assert_msg(max_error<1.e-7,"Largest approximation error %g",max_error);
}
END_TEST

Suite* barriers_suite(void)
{

    Suite* suite = suite_create("Barriers");
    TCase *tc = tcase_create("Positive Orthant");
    tcase_add_test(tc,s_test_pos_orthant);
    tcase_add_test(tc,s_test_pos_orthant_grad);
    tcase_add_test(tc,s_test_pos_orthant_hess);
    tcase_add_test(tc,s_test_pos_orthant_val);
    tcase_add_test(tc,s_test_pos_orthant_nnz);
    tcase_add_test(tc,s_test_pos_orthant_complexity);

    TCase *tce = tcase_create("Exp Cone");
    tcase_add_test(tce,test_exp_nnz);
    tcase_add_test(tce,test_exp_complexity);
    tcase_add_test(tce,test_exp_val);
    tcase_add_test(tce,test_exp_grad);
    tcase_add_test(tce,test_exp_hess);

    suite_add_tcase(suite,tc);
    suite_add_tcase(suite,tce);
    return suite;
}
