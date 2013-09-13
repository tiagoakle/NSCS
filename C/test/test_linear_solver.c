#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "test_linear_solver.h"
#include "linear_solvers.h"
#include "nscs_sp.h"
#include "OpenBLAS/cblas.h"
#include "test_util.h"
#include "smatvec.h"
#include "umfpack.h"

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


    double* x_sol  = (double*)calloc(5,sizeof(double));
    double* b_prod = (double*)calloc(5,sizeof(double));   
    double* x      = (double*)calloc(5,sizeof(double)); 
    
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

START_TEST (test_form_kkt_system)
{
    //printf("Test form KKT system ------------- \n");
     //read H.csv
    int nnzH = 0;
    int n;
    int m;
    read_csv_size("./test/test_data/minentropy_H.csv",&m,&n,&nnzH);
    int* hI= (int*)calloc(nnzH,sizeof(int));
    int* hJ= (int*)calloc(nnzH,sizeof(int));
    double*hV = (double*)calloc(nnzH,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",hI,hJ,hV);   
   // printf("Loaded matrix H of size %i,%i with nnz: %i \n",m,n,nnzH);

    //read A.csv
    int nnzA = 0;  
    read_csv_size("./test/test_data/minentropy_A.csv",&m,&n,&nnzA);
    int* aI= (int*)calloc(nnzA,sizeof(int)); 
    int* aJ= (int*)calloc(nnzA,sizeof(int));

    double*aV = (double*)calloc(nnzA,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",aI,aJ,aV);   
   // printf("Loaded matrix A of size %i,%i with nnz: %i \n",m,n,nnzA);
    
    double mu    = 1.;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    
    int nnzK = m + 2*nnzA+nnzH;

     csi* Ai = NULL;
     csi* Ap = NULL;
     double* Av = NULL;

     int res; 
     res =  form_kkt_system(mu, hI, hJ,hV,nnzH,aI,aJ,aV,
                                     nnzA,\
                                     m,\
                                     n,\
                                     delta,\
                                     gamma,\
                                     &Ai,\
                                     &Ap,\
                                     &Av);

     //printf("Form KKT SYTEM %i\n",res);
     
     int* rI = (int*)calloc(nnzK,sizeof(int));
     int* rJ = (int*)calloc(nnzK,sizeof(int));
     double* rV = (double*)calloc(nnzK,sizeof(double));
     read_csv_triplets("./test/test_data/matlab_built_K.csv",rI,rJ,rV);   
    // printf("Loaded matrix K of size %i,%i with nnz: %i \n",m+n,n+n,nnzK);

     csi* rAi = (csi*)calloc(nnzK,sizeof(csi));
     csi* rAp = (csi*)calloc(n+m+1,sizeof(csi));
     double* rAv = (double*)calloc(nnzK,sizeof(double));
     //Build the csr form of the matlab system
     res  =  umfpack_di_triplet_to_col(n+m,n+m,nnzK,rI,rJ,rV,rAp,rAi,rAv,NULL);
    // printf("Umfpack triplet to col %i \n",res);

     //Compare both systems
      double t1,t2;
      double a1  = 0;
      double a2  = 0;
      double a3  = 0;
    int i = 0; 
     for(i=0;i<nnzK;i++)
     {
        t1 = rAi[i]-Ai[i];    
        t2 = rAv[i]-Av[i];
        a1 += t1*t1;
        a2 += t2*t2;
     }
     for(i=0;i<n+m+1;i++)
     {
        t1 = rAp[i]-Ap[i];
        a3 += t1*t1; 
     }
     //printf("Differences %lf, %lf, %lf\n",a1,a2,a3);
     //fflush(stdout);
     ck_assert_msg(a1<1.e-7,"%g\n",a1);
     ck_assert_msg(a2<1.e-7,"%g\n",a2);
     ck_assert_msg(a3<1.e-7,"%g\n",a3);
}
END_TEST;

//XXX:Something wrong with check hangs this test...
START_TEST(test_load_csv_and_solve_system)
{
    //read H.csv
    int nnzH = 0;
    int n;
    int m;
    read_csv_size("./test/test_data/minentropy_H.csv",&m,&n,&nnzH);
    int* hI= (int*)calloc(nnzH,sizeof(int));
    int* hJ= (int*)calloc(nnzH,sizeof(int));
    double*hV = (double*)calloc(nnzH,sizeof(double));

    read_csv_triplets("./test/test_data/minentropy_H.csv",hI,hJ,hV);   
    //printf("Loaded matrix H of size %i,%i with nnz: %i \n",m,n,nnzH);

    //read A.csv
    int nnzA = 0;  
    read_csv_size("./test/test_data/minentropy_A.csv",&m,&n,&nnzA);
    int* aI= (int*)calloc(nnzA,sizeof(int));
    int* aJ= (int*)calloc(nnzA,sizeof(int));
    double*aV = (double*)calloc(nnzA,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",aI,aJ,aV);   
    //printf("Loaded matrix A of size %i,%i with nnz: %i \n",m,n,nnzA);

    //read rhs.csv 
    int rhs_size;
    rhs_size = read_size("./test/test_data/minentropy_rhs.csv");
    double* rhs = (double*)calloc(rhs_size,sizeof(double));
    read_vector("./test/test_data/minentropy_rhs.csv",rhs);
    //printf("Loaded vector rhs of size %i \n",rhs_size);
    
    //Now call the solver
    double mu = 1;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    
    double* sol = (double*)calloc(n+m,sizeof(double));
//    build_and_solve_linear_system(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, sol, rhs);

     int* Ai = NULL;
     int* Ap = NULL;
     double* Av = NULL;
     void* Numeric = NULL;
 
      int ret = form_kkt_system(mu, hI, hJ, hV,nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, &Ai, &Ap, &Av);
     // printf("Form KKT SYSTEM %i\n",ret);     
      if (ret != 0)
         sol[0] = -1000+ret;     
     
     ret = factor_kkt_system(&Numeric, Ai, Ap, Av, n, m);
    // printf("Factor KKT SYSTEM %i\n",ret);     
    //fflush(stdout);
      if (ret != 0)
         sol[1] = -2000+ret;
      ret = solve_factored_system(Numeric, Ai, Ap, Av, rhs, sol);
    // printf("Solve KKT SYSTEM %i\n",ret);     
      if (ret != 0)
         sol[2] = -3000+ret;

    //Save the solution 
    // write_vector_to_csv("./test/test_data/test_sol.csv",sol,n+m);
    //printf("Saving calculated solution of size %i to csv\n",n+m); 
    //read solution.csv
    int nnz = 0;
    nnz = read_size("./test/test_data/minentropy_sol.csv");
    double* somat = (double*)calloc(nnz,sizeof(double));
    read_vector("./test/test_data/minentropy_sol.csv",somat); 
    //printf("Loaded stored solution of size %i \n",nnz);
    
    double naccum = 0;
    int i = 0;
    double t = 0;
    for(i=0;i<nnz;i++)
    {
        t = (somat[i]-sol[i]);
        naccum += (t*t);
    }
    naccum = sqrt(naccum);

    //since we only keep few digits of accuracy storing the solution in a text file, we cant do better
    ck_assert_msg(naccum<1.e-5,"Norm of difference between calculated and stored sol %g\n",naccum);
    //printf("Difference between saved and calculated sol %lf\n",naccum);    
    //fflush(stdout);

}
END_TEST

//Test the solution to the full HSD system
START_TEST(test_full_hsd_soltuion)
{
    //Reads the system stored in the files
    //minentropy_H,minentropy_A,minentropy_b,..., etc 
    //builds the HSD system and solves it 
    //saves the solution in test_data/debug_Atdy_,...,test_data/debug_dx_solution... etc.
    //Computes the residual and calculates the residual norm for each equation.

    //printf("---------Test of full solve-------\n");

    //read H.csv
    spmat H; 
    read_csv_size("./test/test_data/minentropy_H.csv",&H.m,&H.n,&H.nnz);
    H.I= (int*)calloc(H.nnz,sizeof(int));
    H.J= (int*)calloc(H.nnz,sizeof(int));
    H.V = (double*)calloc(H.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",H.I,H.J,H.V);   
    //printf("Loaded matrix H of size %i,%i with nnz: %i \n",H.m,H.n,H.nnz);

    //read A.csv
    spmat A;
    read_csv_size("./test/test_data/minentropy_A.csv",&A.m,&A.n,&A.nnz);
    A.I = (int*)calloc(A.nnz,sizeof(int));
    A.J = (int*)calloc(A.nnz,sizeof(int));
    A.V = (double*)calloc(A.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",A.I,A.J,A.V);   
    //printf("Loaded matrix A of size %i,%i with nnz: %i \n",A.m,A.n,A.nnz);

    vec b;
    vec c;
    
    char* nameb="./test/test_data/minentropy_b.csv";
    b.n = read_size(nameb);
    b.v = (double*)calloc(b.n,sizeof(double));
    read_vector(nameb,b.v); 
    //printf("Loaded vector b of size %i \n",b.n);

    char* namec="./test/test_data/minentropy_c.csv";
    c.n = read_size(namec);
    c.v = (double*)calloc(c.n,sizeof(double));
    read_vector(namec,c.v); 
    //printf("Loaded vector c of size %i \n",c.n);
    
    //Read the vector that contains [mu,tau,kappa,r3,r4]
    char* name = "./test/test_data/minentropy_doubles.csv";
    vec doubles;
    doubles.n = 5;
    doubles.v = (double*)calloc(5,sizeof(double));
    read_vector(name,doubles.v);
  
    //Allocate space for the residuals and read
    vec *r1 = calloc_vec(A.m);
    name = "./test/test_data/minentropy_r1.csv";
    read_vector(name,r1->v);

    vec *r2 = calloc_vec(A.n);
    name = "./test/test_data/minentropy_r2.csv";
    read_vector(name,r2->v);
    
    //The vector called r5 in the matlab version 
    //Corresponds to the vector r4 here
    vec *r4 = calloc_vec(H.n);
    name = "./test/test_data/minentropy_r5.csv";
    read_vector(name,r4->v);
   
    //Allocate space for the return values
    vec *dx = calloc_vec(A.n);
    vec *dy = calloc_vec(A.m);
    vec *ds = calloc_vec(H.n);
    double dt, dk;
    
    double mu  = doubles.v[0];
    double tau = doubles.v[1];
    double kappa = doubles.v[2];
    double r3    = doubles.v[3];
    double r5    = doubles.v[4];

    double delta = 1.e-10;
    double gamma = 1.e-10;
    int res;
    
    //printf("Loaded residuals ready to call solver\n");

    res  = solve_kkt_system(mu,\
                     H,\
                     A,\
                     b,\
                     c,\
                     tau,\
                     kappa,\
                     delta,\
                     gamma,\
                     *r1,\
                     *r2,\
                     r3,\
                     *r4,\
                     r5,\
                     *dy,\
                     *dx,\
                     &dt,\
                     *ds,\
                     &dk );
     
     //Generate a CRS version of H 
    csi* Hi =(csi*)calloc(H.nnz,sizeof(csi));
    csi* Hp =(csi*)calloc(H.n+1,sizeof(csi)); //XXX: With free variables this can be smaller
    double* Hv =(double*)calloc(H.nnz,sizeof(double));  
    res = umfpack_di_triplet_to_col(H.m,H.n,H.nnz,H.I,H.J,H.V,Hp,Hi,Hv,NULL);     
    //printf("H to CRS %i\n",res); 
 
    //Generate a CRS version of A 
    csi* Ai =(csi*)calloc(A.nnz,sizeof(csi));
    csi* Ap =(csi*)calloc(A.n+1,sizeof(csi)); //XXX: With free variables this can be smaller
    double* Av =(double*)calloc(A.nnz,sizeof(double));  
    res = umfpack_di_triplet_to_col(A.m,A.n,A.nnz,A.I,A.J,A.V,Ap,Ai,Av,NULL);      
    //printf("A to CRS %i\n",res); 

    //Calculate the first residual
    double* res_1 = (double*)calloc(A.m,sizeof(double)); 
    dspmv(A.m,A.n,1.,Ap,Ai,Av,res_1,dx->v);
    cblas_daxpy(A.m,-dt,b.v,1,res_1,1);
    cblas_daxpy(A.m,-1,r1->v,1,res_1,1);
    double n_res_1 = cblas_ddot(A.m,res_1,1,res_1,1); 
    //printf("Norm residual 1 %lf\n",n_res_1) ;
    ck_assert_msg(n_res_1<1e-14,"norm squared of residual 1 exceeded 1e-14; %g",n_res_1);
    
    //Calculate the second residual 
    double* res_2 = (double*)calloc(A.n,sizeof(double));
    dsTpmv(A.m,A.n,-1.,Ap,Ai,Av,res_2,dy->v);
    cblas_daxpy(A.n,dt,c.v,1,res_2,1);
    cblas_daxpy(A.n,-1.,ds->v,1,res_2,1);
    cblas_daxpy(A.n,-1.,r2->v,1,res_2,1);
    double n_res_2 = cblas_ddot(A.n,res_2,1,res_2,1);
    //printf("Norm residual 2 %lf \n",n_res_2);
    ck_assert_msg(n_res_2<1e-14,"norm squared of residual 2 exceeded 1e-14; %g",n_res_2);

    //Calculate the 3rd residual
    double n_res_3 ;
    n_res_3 = cblas_ddot(A.m,b.v,1,dy->v,1) - cblas_ddot(A.n,c.v,1,dx->v,1) -dk -r3;
    //printf("Norm residual 3 %lf \n",n_res_3);
    ck_assert_msg(n_res_3<1e-13,"norm squared of residual 3 exceeded 1e-13; %g",n_res_3);

    //Calculate the 4th residual
    double* res_4 = (double*)calloc(A.n,sizeof(double));
    dspmv(H.m,H.n,mu,Hp,Hi,Hv,res_4,dx->v);
    cblas_daxpy(H.m,1.,ds->v,1,res_4,1);
    cblas_daxpy(H.m,-1.,r4->v,1,res_4,1);
    double n_res_4 = cblas_ddot(H.n,res_4,1,res_4,1);
    //printf("Norm residual 4 %lf\n",n_res_4);
    ck_assert_msg(n_res_4<1e-13,"norm squared of residual 4 exceeded 1e-13; %g",n_res_4);
   
    //Calculate the 5th residual
    double n_res_5 = kappa*dt+tau*dk-r5;
    ck_assert_msg(n_res_5<1e-13,"norm squared of residual 5 exceeded 1e-13; %g",n_res_5);
    //printf("Norm residual 5 %lf \n",n_res_5);

    //Debug stuff
    //write_vector_to_csv("./test/test_data/debug_dy_solution.csv",dy->v,dy->n);
    //write_vector_to_csv("./test/test_data/debug_Atdy_solution.csv",res_2,A.n); 
    //write_vector_to_csv("./test/test_data/debug_ds_solution.csv",ds->v,A.n);
    //write_vector_to_csv("./test/test_data/debug_dt_solution.csv",&dt,1);
    //write_vector_to_csv("./test/test_data/debug_dx_solution.csv",dx->v,dx->n);
    //write_vector_to_csv("./test/test_data/debug_dk_solution.csv",&dk,1);
    //write_vector_to_csv("./test/test_data/debug_r4_solution.csv",r4->v,r4->n);
       
   
}
END_TEST

//Test the solution to the full HSD system using the call 
//that requires no structs
START_TEST(test_full_hsd_soltuion_no_structs)
{
    //Reads the system stored in the files
    //minentropy_H,minentropy_A,minentropy_b,..., etc 
    //builds the HSD system and solves it 
    //Computes the residual and calculates the residual norm for each equation.

    //printf("---------Test of full solve no structs------\n");
    //read H.csv
    spmat H; 
    read_csv_size("./test/test_data/minentropy_H.csv",&H.m,&H.n,&H.nnz);
    H.I= (int*)calloc(H.nnz,sizeof(int));
    H.J= (int*)calloc(H.nnz,sizeof(int));
    H.V = (double*)calloc(H.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_H.csv",H.I,H.J,H.V);   
    //printf("Loaded matrix H of size %i,%i with nnz: %i \n",H.m,H.n,H.nnz);

    //read A.csv
    spmat A;
    read_csv_size("./test/test_data/minentropy_A.csv",&A.m,&A.n,&A.nnz);
    A.I = (int*)calloc(A.nnz,sizeof(int));
    A.J = (int*)calloc(A.nnz,sizeof(int));
    A.V = (double*)calloc(A.nnz,sizeof(double));
    read_csv_triplets("./test/test_data/minentropy_A.csv",A.I,A.J,A.V);   
    //printf("Loaded matrix A of size %i,%i with nnz: %i \n",A.m,A.n,A.nnz);

    vec b;
    vec c;
    
    char* nameb="./test/test_data/minentropy_b.csv";
    b.n = read_size(nameb);
    b.v = (double*)calloc(b.n,sizeof(double));
    read_vector(nameb,b.v); 
    //printf("Loaded vector b of size %i \n",b.n);

    char* namec="./test/test_data/minentropy_c.csv";
    c.n = read_size(namec);
    c.v = (double*)calloc(c.n,sizeof(double));
    read_vector(namec,c.v); 
    //printf("Loaded vector c of size %i \n",c.n);
    
    //Read the vector that contains [mu,tau,kappa,r3,r4]
    char* name = "./test/test_data/minentropy_doubles.csv";
    vec doubles;
    doubles.n = 5;
    doubles.v = (double*)calloc(5,sizeof(double));
    read_vector(name,doubles.v);
  
    //Allocate space for the residuals and read
    vec *r1 = calloc_vec(A.m);
    name = "./test/test_data/minentropy_r1.csv";
    read_vector(name,r1->v);

    vec *r2 = calloc_vec(A.n);
    name = "./test/test_data/minentropy_r2.csv";
    read_vector(name,r2->v);
    
    //The vector called r5 in the matlab version 
    //Corresponds to the vector r4 here
    vec *r4 = calloc_vec(H.n);
    name = "./test/test_data/minentropy_r5.csv";
    read_vector(name,r4->v);
   
    //Allocate space for the return values
    vec *dx = calloc_vec(A.n);
    vec *dy = calloc_vec(A.m);
    vec *ds = calloc_vec(H.n);
    double dt, dk;
    
    double mu  = doubles.v[0];
    double tau = doubles.v[1];
    double kappa = doubles.v[2];
    double r3    = doubles.v[3];
    double r5    = doubles.v[4];

    double delta = 1.e-10;
    double gamma = 1.e-10;
    int res;
    
    //printf("Loaded residuals ready to call solver\n");
    int m;
    int n;
    m   = A.m;
    n   = A.n;

 solve_kkt_system_no_structs(m,n,
                             mu,\
                             H.I,\
                             H.J,\
                             H.V,\
                             H.nnz,\
                             A.I,\
                             A.J,\
                             A.V,\
                             A.nnz,\
                             b.v,\
                             c.v,\
                             tau,\
                             kappa,\
                             delta,\
                             gamma,\
                             r1->v,\
                             r2->v,\
                             r3,\
                             r4->v,\
                             r5,\
                             dy->v,\
                             dx->v,\
                             &dt,\
                             ds->v,\
                             &dk );
    //Generate a CRS version of H 
    csi* Hi =(csi*)calloc(H.nnz,sizeof(csi));
    csi* Hp =(csi*)calloc(H.n+1,sizeof(csi)); //XXX: With free variables this can be smaller
    double* Hv =(double*)calloc(H.nnz,sizeof(double));  
    res = umfpack_di_triplet_to_col(H.m,H.n,H.nnz,H.I,H.J,H.V,Hp,Hi,Hv,NULL);     
    //printf("H to CRS %i\n",res); 
 
    //Generate a CRS version of A 
    csi* Ai =(csi*)calloc(A.nnz,sizeof(csi));
    csi* Ap =(csi*)calloc(A.n+1,sizeof(csi)); //XXX: With free variables this can be smaller
    double* Av =(double*)calloc(A.nnz,sizeof(double));  
    res = umfpack_di_triplet_to_col(A.m,A.n,A.nnz,A.I,A.J,A.V,Ap,Ai,Av,NULL);      
    //printf("A to CRS %i\n",res); 

    //Calculate the first residual
    double* res_1 = (double*)calloc(A.m,sizeof(double)); 
    dspmv(A.m,A.n,1.,Ap,Ai,Av,res_1,dx->v);
    cblas_daxpy(A.m,-dt,b.v,1,res_1,1);
    cblas_daxpy(A.m,-1,r1->v,1,res_1,1);
    double n_res_1 = cblas_ddot(A.m,res_1,1,res_1,1); 
    //printf("Norm residual 1 %lf\n",n_res_1) ;
    ck_assert_msg(n_res_1<1e-14,"norm squared of residual 1 exceeded 1e-14; %g",n_res_1);
    
    //Calculate the second residual 
    double* res_2 = (double*)calloc(A.n,sizeof(double));
    dsTpmv(A.m,A.n,-1.,Ap,Ai,Av,res_2,dy->v);
    cblas_daxpy(A.n,dt,c.v,1,res_2,1);
    cblas_daxpy(A.n,-1.,ds->v,1,res_2,1);
    cblas_daxpy(A.n,-1.,r2->v,1,res_2,1);
    double n_res_2 = cblas_ddot(A.n,res_2,1,res_2,1);
    //printf("Norm residual 2 %lf \n",n_res_2);
    ck_assert_msg(n_res_2<1e-14,"norm squared of residual 2 exceeded 1e-14; %g",n_res_2);

    //Calculate the 3rd residual
    double n_res_3 ;
    n_res_3 = cblas_ddot(A.m,b.v,1,dy->v,1) - cblas_ddot(A.n,c.v,1,dx->v,1) -dk -r3;
    //printf("Norm residual 3 %lf \n",n_res_3);
    ck_assert_msg(n_res_3<1e-13,"norm squared of residual 3 exceeded 1e-13; %g",n_res_3);

    //Calculate the 4th residual
    double* res_4 = (double*)calloc(A.n,sizeof(double));
    dspmv(H.m,H.n,mu,Hp,Hi,Hv,res_4,dx->v);
    cblas_daxpy(H.m,1.,ds->v,1,res_4,1);
    cblas_daxpy(H.m,-1.,r4->v,1,res_4,1);
    double n_res_4 = cblas_ddot(H.n,res_4,1,res_4,1);
    //printf("Norm residual 4 %lf\n",n_res_4);
    ck_assert_msg(n_res_4<1e-13,"norm squared of residual 4 exceeded 1e-13; %g",n_res_4);
   
    //Calculate the 5th residual
    double n_res_5 = kappa*dt+tau*dk-r5;
    ck_assert_msg(n_res_5<1e-13,"norm squared of residual 5 exceeded 1e-13; %g",n_res_5);
    //printf("Norm residual 5 %lf \n",n_res_5);

    //Debug stuff
    //write_vector_to_csv("./test/test_data/debug_dy_solution_no_structs.csv",dy->v,dy->n);
    //write_vector_to_csv("./test/test_data/debug_Atdy_solution_no_structs.csv",res_2,A.n); 
    //write_vector_to_csv("./test/test_data/debug_ds_solution_no_structs.csv",ds->v,A.n);
    //write_vector_to_csv("./test/test_data/debug_dt_solution_no_structs.csv",&dt,1);
    //write_vector_to_csv("./test/test_data/debug_dx_solution_no_structs.csv",dx->v,dx->n);
    //write_vector_to_csv("./test/test_data/debug_dk_solution_no_structs.csv",&dk,1);
    //write_vector_to_csv("./test/test_data/debug_r4_solution_no_structs.csv",r4->v,r4->n);
    //write_vector_to_csv("./test/test_data/debug_r2_solution_no_structs.csv",r2->v,r2->n);
    
}
END_TEST


Suite* linear_solver_suite(void)
{
    Suite* suite = suite_create("Linear system");
    TCase *tc = tcase_create("Case1");
    tcase_set_timeout(tc,4);
    tcase_add_test(tc,test_solve_linear_system);
    tcase_add_test(tc,test_solve_random_system);
    tcase_add_test(tc,test_form_kkt_system);
    tcase_add_test(tc,test_load_csv_and_solve_system);
    tcase_add_test(tc,test_full_hsd_soltuion_no_structs);
    tcase_add_test(tc,test_full_hsd_soltuion);
    suite_add_tcase(suite,tc);
    return suite;
}


