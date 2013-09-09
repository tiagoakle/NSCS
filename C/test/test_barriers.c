//Test the functions that evaluate the barriers
#include "barriers.h"
#include "check.h"
#include "stdio.h"
#include "test_barriers.h"

void test_pos_orthant_feas(void)
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

START_TEST(s_test_pos_orthant)
{
    test_pos_orthant_feas();
}
END_TEST

//Test gradient for the positive orthant barrier 
void test_pos_orthant_grad(void)
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

START_TEST(s_test_pos_orthant_grad)
{
    test_pos_orthant_grad();
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
    suite_add_tcase(suite,tc);
    return suite;
}
