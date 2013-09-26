//Test the functions that evaluate the barriers
#include "barriers.h"
#include "check.h"
#include "stdio.h"
#include "OpenBLAS/cblas.h"
#include "test_barriers.h"
#include "umfpack.h"
#include "smatvec.h"
#include "linear_solvers.h"
#undef I

START_TEST(s_test_pos_orthant)
{
    int n = 100;
    double* x = (double*)calloc(n,sizeof(double)); 

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
    double* x = (double*)calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1+i*1.0;
    }
    double* g = (double*)calloc(n,sizeof(double));
    pos_orthant_grad(g,x,n);
    
    for(i=0;i<n;i++)
    {
        if((g[i] - 1./x[i])>1.e-16) break;
    }
    
    ck_assert(i==n); 
    free(x);
    free(g);
}
END_TEST


//Test the hessian for the positive orthant barrier
START_TEST(s_test_pos_orthant_hess)
{
    int n = 100;
    double* x = (double*)calloc(n,sizeof(double)); 

    int i = 0;
    for(i=0;i<n;i++)
    { 
        x[i] = 1+i*1.0;
    }

    double delta = 1;

    csi* HI = (csi*)calloc(n,sizeof(csi));
    csi* HJ = (csi*)calloc(n,sizeof(csi));
    double* HV = (double*)calloc(n,sizeof(double));
    pos_orthant_hessian(HI,HJ,HV,0,x,n,delta);
    
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
    double* x = (double*)calloc(n,sizeof(double)); 

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

START_TEST(test_pos_orthant_feasible)
{
    csi n = 100;
    csi i;
    double* x = (double*)calloc(n,sizeof(double));
    for(i=0;i<n;i++) x[i] = 1.0;
    ck_assert(pos_orthant_feas(x,n));
    x[n-1]=-1;
    ck_assert(!pos_orthant_feas(x,n));
}
END_TEST

START_TEST(test_exp_primal_feasible)
{
    csi n = 100;
    csi i,j;
    int tests = 1000;
    double x[3];
    for(j=0;j<tests;j++)
    { 
        x[0] = rand()/(double)RAND_MAX - 0.5;
        x[2] = rand()/(double)RAND_MAX+1.;
        x[1] = x[2]*exp(x[0]/x[2]) + rand()/(double)RAND_MAX;
        ck_assert_msg(exp_primal_feas(x),"Test retured infeasible x2e(x0/x2):%g, <x1: %g, (x1:%g,x2:%g,x3:%g)",-x[2]*exp(x[0]/x[2]),x[1],x[0],x[1],x[2]);

        x[0] = rand()/(double)RAND_MAX - 0.5;
        x[2] = rand()/(double)RAND_MAX + 1.;
        x[1] = x[2]*exp(x[0]/x[2]) - rand()/(double)RAND_MAX;

        ck_assert_msg(!exp_primal_feas(x),"Test retured feasible x2e(x0/x2):%g, >x1: %g, (x1:%g,x2:%g,x3:%g)",-x[2]*exp(x[0]/x[2]),x[1],x[0],x[1],x[2]);

    }
}
END_TEST


START_TEST(test_exp_dual_feasible)
{
    csi n = 100;
    csi i,j;
    int tests = 1000;
    double x[3];
    for(j=0;j<tests;j++)
    { 
        x[0] = -rand()/(double)RAND_MAX-1;
        x[2] = rand()/(double)RAND_MAX -0.5;
        x[1] = -x[0]*exp(x[2]/x[0]-1.) + rand()/(double)RAND_MAX;
        ck_assert_msg(exp_dual_feas(x),"Point %g,%g,%g reported dual infeasible -x[0]*exp(x[2]/x[0]-1.) %g, x[2] %g ",x[0],x[1],x[2],-x[0]*exp(x[2]/x[0]-1.),x[1]);
        
        x[0] = -rand()/(double)RAND_MAX-1;
        x[2] = rand()/(double)RAND_MAX -0.5;
        x[1] = -x[0]*exp(x[2]/x[0]-1.) - rand()/(double)RAND_MAX;
        ck_assert_msg(!exp_dual_feas(x),"Point %g,%g,%g reported dual feasible -x[0]*exp(x[2]/x[0]-1.) %g, x[2] %g ",x[0],x[1],x[2],-x[0]*exp(x[2]/x[0]-1.),x[1]);
    }
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


START_TEST(test_exp_hess_numerical)
{
   csi HI[9];
   csi HJ[9];
   double HV[9];
   double x[] = {0.5,exp(1),1.};
   double xt[3];
   double g[3];
   double gt[3];
   double Hs[3];
   exp_hessian(HI,HJ,HV,0,x,0.); 

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
    int tests = 100;
    int t = 0;
    int i = 0;
    double a = 1; 
    double s[3];
    double n_sq;
    double s_norm;

    double max_error = 0.;
    double min_err_ap;
    //evaluate the gradient
    exp_grad(g,x);
    for(t=0;t<tests;t++)
    {
        

        //Build a random direction 
        for(i=0;i<3;i++)
            s[i] = rand() / (double)RAND_MAX - 0.5;
        //Normalize
        s_norm = cblas_ddot(3,s,1,s,1);
        s_norm = sqrt(s_norm);
        cblas_dscal(3,1./s_norm,s,1);
        //Evaluate the product Hs
        dspmv(3,3,1.,Ap,Ai,Av,Hs,s); 
        
        //Calculate the approximation
        min_err_ap = 1.; 


        for(a=1.0;a>1e-10;a=a/2){

            cblas_dcopy(3,x,1,xt,1); //copy x into xt
            cblas_daxpy(3,a,s,1,xt,1); //Shift xt by a*s

            //Evaluate the gradients
            exp_grad(gt,xt);
       
            //Test
            Hs[0]  = 0;
            Hs[1]  = 0;
            Hs[2]  = 0;

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
    ck_assert_msg(max_error<1.e-6,"Largest approximation error %g",max_error);
}
END_TEST


START_TEST(test_gradient_values_two_cone)
{
    //Chooses a feasible point and generates the gradient for a two cone problem
    //Then calls the pos_orthant_grad and exp_grad and compares the
    //values of the non zeros
    
    int n = 13; //10 positive and 30 exponential cones
    int k_count = 2;
    csi nK[k_count];
    int tK[k_count];
    csi i;

    //Set the number of variables in each cone
    //and the type of each cone
    nK[0] = 10;
    tK[0] = 0;
    nK[1] = 3;
    tK[1] = 3;

    //Build the problem description
    spmat A; 
    A.n   = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;
    problem.delta = 0.;

    double* g=(double*)calloc(n,sizeof(double)); 
    double* cg=(double*)calloc(n,sizeof(double)); 
      
    //Choose a strictly feasible point
    double* x = (double*)calloc(n,sizeof(double));    
    for(i=0;i<10;i++) x[i] = 0.5;
    //Choose feasible points for the exponential cone
    x[10] = rand()/(double)RAND_MAX - 0.5 ; 
    x[12] = rand()/(double)RAND_MAX + 0.01;
    x[11] = x[12]*exp(x[10]/x[12])+1;
    
   //Evaluate the hessian    
    eval_grad(problem,x,g);
   //Eval the pos orthant hessian
   
   pos_orthant_grad(cg,x,nK[0]);
   exp_grad(cg+nK[0],x+nK[0]);
 
   for(i=0;i<n;i++) ck_assert_msg(g[i]==cg[i],"Gradient g[%i] wrong, has %g should be %g\n",i,g[i],cg[i]);
   free(x);
   free(g);
   free(cg);

}
END_TEST 


START_TEST(test_hessian_assemble)
{ 
    //Generate a hessian with two cones, positive orthant and exp cone
    
    int n = 5; //2 positive vars and 1 exp cone
    int k_count = 2;
    csi nK[k_count];
    int tK[k_count];
    csi i;
    //Set the number of variables in each cone
    //and the type of each cone
    nK[0] = 2;
    tK[0] = 0;
    for(i=1;i<k_count;i++){nK[i] = 3; tK[i]=3;}
    //Build the problem description
    spmat A;
    A.n     = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;
    //Calculate the number of non zeros of the hessian
    int nnzH = hessian_nnz(problem);
    //Allocate space for the hessian
    csi* HI=(csi*)calloc(nnzH,sizeof(csi));
    csi* HJ=(csi*)calloc(nnzH,sizeof(csi));
    double* HV=(double*)calloc(nnzH,sizeof(double)); 
    
    state_t state;
    spmat H; 
    H.I = HI;
    H.J = HJ;
    H.V = HV;
    H.n = n;
    H.m = n;
    H.nnz = nnzH;
    state.H = H;

    //Choose a strictly feasible point
    double* x = (double*)calloc(n,sizeof(double));    
    for(i=0;i<2;i++) x[i] = 0.5;
    //Choose feasible points for the exponential cone
    for(i=2;i<n;i=i+3){ x[i] = rand()/(double)RAND_MAX - 0.5 ; x[i+2] = rand()/(double)RAND_MAX ; x[i+1] = x[i+2]*exp(x[i]/x[i+2])+10;}
    
    state.x = x; 
   //Evaluate the hessian    
    eval_hess(problem,x,state);
   //Compress to csr
   //convert the hessian into compressed column form
   csi*   Ai = (csi*)calloc(nnzH,sizeof(csi));
   csi*   Ap = (csi*)calloc(n+1,sizeof(csi));     
   double* Av = (double*)calloc(nnzH,sizeof(double)); 
   csi* Map=NULL;
   //Execute the call
    int status;
    status  =  umfpack_di_triplet_to_col(n,n,nnzH,HI,HJ,HV,Ap,Ai,Av,Map);
 
   //The resulting indices should be these
   csi rAi[11] = {0,1,2,3,4,2,3,4,2,3,4};
   csi rAp[6]  = {0,1,2,5,8,11}; 
   //Compare 
   for(i=0;i<11;i++) ck_assert_msg(Ai[i]==rAi[i],"Hessian Ai[%i] wrong, has %i should be %i\n",i,Ai[i],rAi[i]);
   for(i=0;i<6;i++) ck_assert_msg(Ap[i]==rAp[i],"Hessian Ap[%i] wrong, has %i should be %i\n",i,Ap[i],rAp[i]);

   free(Ai);
   free(Ap);
   free(Av);
   free(HI);
   free(HJ);
   free(x);

}
END_TEST

START_TEST(test_general_nnz)
//void test_general_hessian()
{
    //Chooses a feasible point and generates the hessian for a multiple
    //cone barrier. It then does a numerical approximation of the hessian
    //allong a series of random directions
    
    int n = 100; //10 positive and 30 exponential cones
    int k_count = 31;
    csi nK[k_count];
    int tK[k_count];
    csi i;
    //Set the number of variables in each cone
    //and the type of each cone
    nK[0] = 10;
    tK[0] = 0;
    for(i=1;i<k_count;i++){nK[i] = 3; tK[i]=3;}

    //Build the problem description
    spmat A;
    A.n = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;

    //Calculate the number of non zeros of the hessian
    int nnzH = hessian_nnz(problem);
    ck_assert_msg(nnzH - 10 -9*30 == 0,"nnzH returned wrong count: %i",nnzH);
}
END_TEST

START_TEST(test_values_two_cone_hessian)
{
    //Chooses a feasible point and generates the hessian for a two cone problem
    //Then calls the pos_orthant_hessian and exp_hessian and compares the
    //values of the non zeros
    
    int n = 13; //10 positive and 30 exponential cones
    int k_count = 2;
    csi nK[k_count];
    int tK[k_count];
    csi i;

    //Set the number of variables in each cone
    //and the type of each cone
    nK[0] = 10;
    tK[0] = 0;
    nK[1] = 3;
    tK[1]=3;

    //Build the problem description
    spmat A;
    A.n = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;
    problem.delta = 0.;

    //Calculate the number of non zeros of the hessian
    int nnzH = hessian_nnz(problem);
    //Allocate space for the hessian
    csi* HI=(csi*)calloc(nnzH,sizeof(csi));
    csi* HJ=(csi*)calloc(nnzH,sizeof(csi));
    double* HV=(double*)calloc(nnzH,sizeof(double)); 
    double* cHV=(double*)calloc(nnzH,sizeof(double)); 
      
    state_t state;
    spmat H; 
    H.I = HI;
    H.J = HJ;
    H.V = HV;
    H.n = n;
    H.m = n;
    H.nnz = nnzH;
    state.H = H;

    //Choose a strictly feasible point
    double* x = (double*)calloc(n,sizeof(double));    
    for(i=0;i<10;i++) x[i] = 0.5;
    //Choose feasible points for the exponential cone
    x[10] = rand()/(double)RAND_MAX - 0.5 ; 
    x[12] = rand()/(double)RAND_MAX + 0.01;
    x[11] = x[12]*exp(x[10]/x[12])+1;
    
    state.x = x;
   
   //Evaluate the hessian    
    eval_hess(problem,x,state);
   //Eval the pos orthant hessian
   pos_orthant_hessian(HI,HJ,cHV,0,x,nK[0],0.);
   exp_hessian(HI,HJ,cHV+nK[0],0,x+nK[0],0.);
    
   for(i=0;i<nnzH;i++) ck_assert_msg(cHV[i]==HV[i],"Hessian HV[%i] wrong, has %g should be %g\n",i,HV[i],cHV[i]);
   free(x);
   free(cHV);
   free(HV);
}
END_TEST 

START_TEST(test_general_call_hessian_gradient_positive_orthant)
{
    //Chooses a feasible point and generates the hessian and gradient for a positive orhtant cone
    //but uses the general call and a point for which the values are known
    
    int n = 10; //10 positive variables
    int k_count = 1;

    csi nK[k_count];
    int tK[k_count];
    csi i;

    nK[0] = 10;
    tK[0] = 0;
   
    //Choose a strictly feasible point
    double* x = (double*)calloc(n,sizeof(double));    
    //Allocate space for the gradient
    double* g = (double*)calloc(n,sizeof(double));    
    for(i=0;i<10;i++) x[i] = 1./((double)i+1);

    //Build the problem description
    spmat A;
    A.n = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;
    problem.delta = 0.;

    //Calculate the number of non zeros of the hessian
    int nnzH = hessian_nnz(problem);
    ck_assert_msg(nnzH==n,"nnzH should be %i, id %i",n,nnzH);
    //Allocate space for the hessian
    csi*    HI=(csi*)calloc(nnzH,sizeof(csi));
    csi*    HJ=(csi*)calloc(nnzH,sizeof(csi));
    double* HV=(double*)calloc(nnzH,sizeof(double)); 
      
    state_t state;
    spmat H; 
    H.I = HI;
    H.J = HJ;
    H.V = HV;
    H.n = n;
    H.m = n;
    H.nnz = nnzH;
    state.H = H;
   
    //Evaluate the gradient 
    eval_grad(problem,x,g);
    //Evaluate the hessian    
    eval_hess(problem,x,state);

    for(i=0;i<nnzH;i++) ck_assert_msg(HV[i]==(i+1.)*(i+1),"Diagonal %i of hessian should be %g is %g",i,(1.+i)*(1.+i),HV[i]);
    for(i=0;i<n;i++)    ck_assert_msg(g[i]==-(i+1.),"Entry %i of gradient should be  %g is %g",i,-(1.+i),g[i]);

    free(x);
    free(g);
    free(HI);
    free(HJ);
    free(HV);
}
END_TEST

START_TEST(test_general_hessian_numerical)
{
    //Chooses a feasible point and generates the hessian for a multiple
    //cone barrier. It then does a numerical approximation of the hessian
    //allong a series of random directions
    
    int n = 3; //10 positive and 30 exponential cones
    int k_count = 1;

    csi nK[k_count];
    int tK[k_count];
    csi i;

    //Set the number of variables in each cone
    //and the type of each cone
    nK[0] = 0;
    tK[0] = 0;
    //Assign the type and size to the exponential cones    
    for(i=1;i<k_count;i++){nK[i] = 3; tK[i]=3;}
    
    double* x = (double*)calloc(n,sizeof(double));    
    for(i=0;i<nK[0];i++) x[i] = 0.5;

    //Choose feasible points for the exponential cone
    for(i=nK[0];i<n;i=i+3){ x[i] = rand()/(double)RAND_MAX - 0.5 ; x[i+2] = rand()/(double)RAND_MAX ; x[i+1] = x[i+2]*exp(x[i]/x[i+2])+1;}
  
    //Build the problem description
    spmat A;
    A.n = n;
    problem_t problem;
    problem.A = A;
    problem.k_count = k_count;
    problem.nK = nK;
    problem.tK = tK;
    problem.delta = 0.;

    //Calculate the number of non zeros of the hessian
    int nnzH = hessian_nnz(problem);
    //Allocate space for the hessian
    csi*    HI=(csi*)calloc(nnzH,sizeof(csi));
    csi*    HJ=(csi*)calloc(nnzH,sizeof(csi));
    double* HV=(double*)calloc(nnzH,sizeof(double)); 
      
    state_t state;
    spmat H; 
    H.I = HI;
    H.J = HJ;
    H.V = HV;
    H.n = n;
    H.m = n;
    H.nnz = nnzH;
    state.H = H;
   
    //Evaluate the hessian    
    eval_hess(problem,x,state);
    
   //convert the hessian into compressed column form
   csi*    Ai = (csi*)calloc(nnzH,sizeof(csi));
   csi*    Ap = (csi*)calloc(n+1,sizeof(csi));     
   double* Av = (double*)calloc(nnzH,sizeof(double)); 
   csi* Map=NULL;
   //Execute the call
   int status;
   status  =  umfpack_di_triplet_to_col(n,n,nnzH,HI,HJ,HV,Ap,Ai,Av,Map);
   ck_assert_msg(status==0,"Unable to generate csr form");

   //Each test makes a product with the hessian and compares with the gradient
   int tests = 100;
   int t = 0;
   double a = 1; 
   //Allocate a vector for the direction
   double* s  = (double*)calloc(n,sizeof(double));
   //Allocate a vector for the product
   double* Hs = (double*)calloc(n,sizeof(double));
   double* sol = (double*)calloc(n,sizeof(double));
   //Allocate vectors for the gradients 
   double* g  = (double*)calloc(n,sizeof(double));
   double* gt = (double*)calloc(n,sizeof(double));
   //Allocate a vector for the trail point
   double* xt = (double*)calloc(n,sizeof(double));
   
   //Evaluate the gradient at x
   eval_grad(problem,x,g);

   double n_sq;
   double s_norm;

   double max_error = 0.;
   double min_err_ap;
   for(t=0;t<tests;t++)
    {
       
        //Build a random direction 
        for(i=0;i<n;i++)
            s[i] = rand() / (double)RAND_MAX - 0.5;
        //Normalize
        s_norm = cblas_ddot(n,s,1,s,1);
        s_norm = sqrt(s_norm);
        cblas_dscal(n,1./s_norm,s,1);

        //Initialize the minimum to a large value 
        min_err_ap = 1.e10; 

             
        for(i=0;i<n;i++){ Hs[i] = 0.; sol[i] = 0.;}

        //Evaluate the product H*s
        dspmv(n,n,1.,Ap,Ai,Av,Hs,s); 
        
        //For debug
        int status = solve_linear_system(sol,HI,HJ,HV,nnzH,Hs,n);
        ck_assert_msg(status==0,"status %i",status);
        cblas_daxpy(n,-1.,s,1,sol,1);
        printf("Check product %f\n",sqrt(cblas_ddot(n,sol,1,sol,1)));

        //Calculate the approximation
        for(a=1.0;a>1e-10;a=a/2){

            cblas_dcopy(n,x,1,xt,1); //copy x into xt
            cblas_daxpy(n,a,s,1,xt,1); //Shift xt by a*s

            //Evaluate the gradients
            eval_grad(problem,xt,gt);
          
            //calculate -1/a(a*Hs+gs-gt)
            //as Hs+(-g+gt)(-1/a)
            cblas_daxpy(n,-1.,g,1,gt,1);
            cblas_daxpy(n,-1./a,gt,1,Hs,1);

            //Calculate the norm squared
            n_sq = cblas_ddot(n,Hs,1,Hs,1);
            min_err_ap = fmin(min_err_ap,sqrt(n_sq));//keep the smallst error

        }
        max_error = fmax(min_err_ap,max_error); //Find the largest over all the random tests
    }
    ck_assert_msg(max_error<1.e-6,"Maximum approximation error %g >1.e-6",max_error);

    //Clean up
    free(x);
    free(HI);
    free(HJ);
    free(HV);
    free(Ai);
    free(Ap);
    free(Av);
    free(g);
    free(gt);
    free(xt);
    free(Hs);
    free(sol);
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
    tcase_add_test(tc,test_pos_orthant_feasible);

    TCase *tce = tcase_create("Exp Cone");
    tcase_add_test(tce,test_exp_nnz);
    tcase_add_test(tce,test_exp_complexity);
    tcase_add_test(tce,test_exp_val);
    tcase_add_test(tce,test_exp_grad);
    tcase_add_test(tce,test_exp_hess_numerical);
    tcase_add_test(tce,test_exp_primal_feasible);
    tcase_add_test(tce,test_exp_dual_feasible);

    TCase *tcf = tcase_create("Milti Cone");
    tcase_add_test(tcf, test_gradient_values_two_cone);
    tcase_add_test(tcf, test_general_hessian_numerical);
    tcase_add_test(tcf, test_hessian_assemble);
    tcase_add_test(tcf, test_general_nnz);
    tcase_add_test(tcf, test_values_two_cone_hessian);
    tcase_add_test(tcf, test_general_call_hessian_gradient_positive_orthant);
    suite_add_tcase(suite,tc);
    suite_add_tcase(suite,tce);
    suite_add_tcase(suite,tcf);

    return suite;
}
