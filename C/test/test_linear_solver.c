#include <stdio.h>
#include "linear_solvers.h"
#include "test_util.h"
#include <stdlib.h>
#include <math.h>
#include "test_linear_solver.h"
#include "nscs_sp.h"

START_TEST (test_solve_linear_system)
{

    //Test to call linear_solvers
    int n = 5 ;
    int Ap [ ] = {0, 2, 5, 9, 10, 12} ; 
    int Aj [ ] = { 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4} ;
    int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
    double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
    double b [ ] = {8., 45., -3., 3., 19.} ;
    double x [5] ;
    double x_sol[] = {1.,2.,3.,4.,5.};


    int status = solve_linear_system(x, Ai, Aj, Ax, 12, b , n);

    double t = 0;
    double acc = 0;
    int i = 0;
    for(i=0;i<5;i++)
    {
        t = (x[i]-x_sol[i]);
        acc += t*t;
    }

    ck_assert(acc<1.e-10);
   
}
END_TEST

START_TEST (test_solve_random_system)
{
    
    //Test to call linear_solvers
    int n = 5 ;
    int Ap [ ] = {0, 2, 5, 9, 10, 12} ; 
    int Aj [ ] = { 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4} ;
    int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
    double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;


    double* x_sol  = calloc(5,sizeof(double));
    double* b_prod = calloc(5,sizeof(double));   
    double* x      = calloc(5,sizeof(double)); 
    
    //Build random rhs 
    int i = 0;
    for(i=0;i<5;i++)
    {
        x_sol[i] = ((double)random())/((double)RAND_MAX);
    }
    //Calculate the rhs 
    dspmv(5,5,1.,Ap,Ai,Ax,b_prod,x_sol);
    
    //Solve the system
    int status = solve_linear_system(x, Ai, Aj, Ax, 12, b_prod , n);

    double t = 0;
    double acc = 0;
    for(i=0;i<5;i++)
    {
        t = (x[i]-x_sol[i]);
        acc += t*t;
    }

    ck_assert(acc<1.e-10);
   
}
END_TEST

//void test_load_triplets()
//{
//    //read H.csv
//    int nnz = 0;
//    int n;
//    int m;
//    read_csv_size("./test/H.csv",&n,&m,&nnz);
//    printf("Loaded matrix of size %i,%i with nnz: %i \n",m,n,nnz);
//    
//    int* I= calloc(nnz,sizeof(int));
//    int* J= calloc(nnz,sizeof(int));
//    double*V = calloc(nnz,sizeof(double));
//    read_csv_triplets("./test/H.csv",I,J,V);
//    //Write down a copy
//    write_matrix_to_csv("./test/H2.csv",I,J,V,m,n,nnz); 
//}

START_TEST (form_kkt_system_in_triplet_format)
{
     //read H.csv
    int nnzH = 0;
    int n;
    int m;
    read_csv_size("./test/test_data/minentropy_H.csv",&m,&n,&nnzH);
    int* hI= calloc(nnzH,sizeof(int));
    int* hJ= calloc(nnzH,sizeof(int));
    double*hV = calloc(nnzH,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",hI,hJ,hV);   
    printf("Loaded matrix H of size %i,%i with nnz: %i \n",m,n,nnzH);

    //read A.csv
    int nnzA = 0;  
    read_csv_size("./test/test_data/minentropy_A.csv",&m,&n,&nnzA);
    int* aI= calloc(nnzA,sizeof(int));
    int* aJ= calloc(nnzA,sizeof(int));
    double*aV = calloc(nnzA,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",aI,aJ,aV);   
    printf("Loaded matrix A of size %i,%i with nnz: %i \n",m,n,nnzA);
    
    double mu    = 1.;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    
    int nnzK = m + 2*nnzA+nnzH;

     int* I = calloc(nnzK,sizeof(int));
     int* J = calloc(nnzK,sizeof(int));
     double* V = calloc(nnzK,sizeof(double));
     int res; 
     res =  form_kkt_system_triplets(mu, hI, hJ,hV,nnzH,aI,aJ,aV,
                                     nnzA,\
                                     m,\
                                     n,\
                                     delta,\
                                     gamma,\
                                     I,\
                                     J,\
                                     V);

     
     int* rI = calloc(nnzK,sizeof(int));
     int* rJ = calloc(nnzK,sizeof(int));
     double* rV = calloc(nnzK,sizeof(double));
     read_csv_triplets("./test/test_data/matlab_built_K.csv",rI,rJ,rV);   
     printf("Loaded matrix K of size %i,%i with nnz: %i \n",m+n,n+n,nnzK);

     //Now check that the constructed and the saved K match
    int i;
    double errI =0;
    double errJ =0;
    double errV =0;
    double t;
    for(i=0;i<nnzK;i++)
    { 
       t = I[i] - rI[i];
       errI += t*t; 
       t = J[i] - rJ[i];
       errI += t*t;
       t = V[i] - rV[i];
       errV += t*t;
    }
    printf("Built KKT system vs stored\n");
    printf("Squared error in I %lf, in J %lf, in V %lf\n",errI,errJ,errV); 
}
END_TEST;

//XXX:Something wrong with check hangs this test...
void test_load_csv_and_solve_system(void)
{
    //read H.csv
    int nnzH = 0;
    int n;
    int m;
    read_csv_size("./test/test_data/minentropy_H.csv",&m,&n,&nnzH);
    int* hI= calloc(nnzH,sizeof(int));
    int* hJ= calloc(nnzH,sizeof(int));
    double*hV = calloc(nnzH,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",hI,hJ,hV);   
    printf("Loaded matrix H of size %i,%i with nnz: %i \n",m,n,nnzH);

    //read A.csv
    int nnzA = 0;  
    read_csv_size("./test/test_data/minentropy_A.csv",&m,&n,&nnzA);
    int* aI= calloc(nnzA,sizeof(int));
    int* aJ= calloc(nnzA,sizeof(int));
    double*aV = calloc(nnzA,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",aI,aJ,aV);   
    printf("Loaded matrix A of size %i,%i with nnz: %i \n",m,n,nnzA);

    //read rhs.csv 
    int rhs_size;
    rhs_size = read_size("./test/test_data/minentropy_rhs.csv");
    double* rhs = calloc(rhs_size,sizeof(double));
    read_vector("./test/test_data/minentropy_rhs.csv",rhs);
    printf("Loaded vector rhs of size %i \n",rhs_size);
    
    //Now call the solver
    double mu = 1;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    
    double* sol = calloc(n+m,sizeof(double));
//    build_and_solve_linear_system(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, sol, rhs);

     int* Ai = NULL;
     int* Ap = NULL;
     double* Av = NULL;
     void* Numeric = NULL;
 
      int ret = form_kkt_system(mu, hI, hJ, hV,nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, &Ai, &Ap, &Av);
      printf("Form KKT SYSTEM %i\n",ret);     
      if (ret != 0)
         sol[0] = -1000+ret;     
     
     ret = factor_kkt_system(&Numeric, Ai, Ap, Av, n+m);
     printf("Factor KKT SYSTEM %i\n",ret);     
     fflush(stdout);
      if (ret != 0)
         sol[1] = -2000+ret;
      ret = solve_factored_system(Numeric, Ai, Ap, Av, rhs, sol);
     printf("Solve KKT SYSTEM %i\n",ret);     
      if (ret != 0)
         sol[2] = -3000+ret;

    //Save the solution 
    write_vector_to_csv("./test/test_data/test_sol.csv",sol,n+m);
    printf("Saving calculated solution of size %i to csv\n",n+m); 
    //read solution.csv
    int nnz = 0;
    nnz = read_size("./test/test_data/minentropy_sol.csv");
    double* somat = calloc(nnz,sizeof(double));
    read_vector("./test/test_data/minentropy_sol.csv",somat); 
    printf("Loaded stored solution of size %i \n",nnz);
    
    double naccum = 0;
    int i = 0;
    double t = 0;
    for(i=0;i<nnz;i++)
    {
        t = (somat[i]-sol[i]);
        naccum += (t*t);
    }
    naccum = sqrt(naccum);

    //ck_assert(naccum<1.e-10);
    printf("Difference between saved and calculated sol %lf\n",naccum);    
    fflush(stdout);

}

//Test the solution to the full HSD system
void test_full_hsd_soltuion(void)
{
    //read H.csv
    spmat H; 
    read_csv_size("./test/test_data/minentropy_H.csv",&H.m,&H.n,&H.nnz);
    H.I= calloc(H.nnz,sizeof(int));
    H.J= calloc(H.nnz,sizeof(int));
    H.V = calloc(H.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",H.I,H.J,H.V);   
    printf("Loaded matrix H of size %i,%i with nnz: %i \n",H.m,H.n,H.nnz);

    //read A.csv
    spmat A;
    read_csv_size("./test/test_data/minentropy_A.csv",&A.m,&A.n,&A.nnz);
    A.I = calloc(A.nnz,sizeof(int));
    A.J = calloc(A.nnz,sizeof(int));
    A.V = calloc(A.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",A.I,A.J,A.V);   
    printf("Loaded matrix A of size %i,%i with nnz: %i \n",A.m,A.n,A.nnz);

    vec b;
    vec c;
    
    char* nameb="./test/test_data/minentropy_b.csv";
    b.n = read_size(nameb);
    read_vector(nameb,b.v); 
    printf("Loaded vector b of size %i \n",b.n);

    char* namec="./test/test_data/minentropy_c.csv";
    c.n = read_size(namec);
    read_vector(namec,c.v); 
    printf("Loaded vector c of size %i \n",c.n);
    //This??
    free(nameb);
    free(namec);

    //Allocate space for the residuals
    vec *r1 = calloc_vec(A.m);
    vec *r2 = calloc_vec(A.n);
    double *r3;
    vec *r4 = calloc_vec(H.n);
    double *r5;
   
    //Allocate space for the return values
    vec *dx = calloc_vec(A.n);
    vec *dy = calloc_vec(A.m);
    vec *ds = calloc_vec(H.n);
    double dt, dk;

    double tau = 10;
    double kappa = 10;
    double mu = 1.;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    int res;
    res  = solve_kkt_system(mu,\
                     H,\
                     A,\
                     b,\
                     c,\
                     tau,\
                     kappa,\
                     delta,\
                     gamma,\
                     r1,\
                     r2,\
                     r3,\
                     r4,\
                     r5,\
                     *dy,\
                     *dx,\
                     &dt,\
                     *ds,\
                     &dk );

 
}


Suite* linear_solver_suite(void)
{
    Suite* suite = suite_create("Linear system");
    TCase *tc = tcase_create("Case1");
    tcase_set_timeout(tc,4);
    tcase_add_test(tc,test_solve_linear_system);
    tcase_add_test(tc,test_solve_random_system);
    tcase_add_test(tc,form_kkt_system_in_triplet_format);
//    tcase_add_test(tc,test_load_csv_and_solve_system);
    suite_add_tcase(suite,tc);
    return suite;
}


