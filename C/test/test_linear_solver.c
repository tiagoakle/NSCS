#include <stdio.h>
#include "linear_solvers.h"
#include "test_util.h"
#include <stdlib.h>
#include <math.h>
//Test to call linear_solvers
int n = 5 ;
int Ap [ ] = {0, 2, 5, 9, 10, 12} ; 


int Aj [ ] = { 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4} ;
int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double b [ ] = {8., 45., -3., 3., 19.} ;
double x [5] ;

int test_solve_linear_system()
{
    int status = solve_linear_system(x, Ai, Aj, Ax, 12, b , n);
    size_t i = 0;
    printf("RHS \n");
    for(i=0;i<n;i++)
        printf("%f\n",b[i]);

    printf("Solution \n");
    for(i=0;i<n;i++)
        printf("%f\n",x[i]);
}

void test_load_triplets()
{
    //read H.csv
    int nnz = 0;
    int n;
    int m;
    read_csv_size("./test/H.csv",&n,&m,&nnz);
    printf("Loaded matrix of size %i,%i with nnz: %i \n",m,n,nnz);
    
    int* I= calloc(nnz,sizeof(int));
    int* J= calloc(nnz,sizeof(int));
    double*V = calloc(nnz,sizeof(double));
    read_csv_triplets("./test/H.csv",I,J,V);
    //Write down a copy
    write_matrix_to_csv("./test/H2.csv",I,J,V,m,n,nnz); 
}

void test_load_csv_and_solve_system()
{
    //read H.csv
    int nnzH = 0;
    int n;
    int m;
    read_csv_size("./test/H.csv",&m,&n,&nnzH);
    int* hI= calloc(nnzH,sizeof(int));
    int* hJ= calloc(nnzH,sizeof(int));
    double*hV = calloc(nnzH,sizeof(double));
    read_csv_triplets("./test/H.csv",hI,hJ,hV);   
    printf("Loaded matrix H of size %i,%i with nnz: %i \n",m,n,nnzH);

    //read A.csv
    int nnzA = 0;  
    read_csv_size("./test/A.csv",&m,&n,&nnzA);
    int* aI= calloc(nnzA,sizeof(int));
    int* aJ= calloc(nnzA,sizeof(int));
    double*aV = calloc(nnzA,sizeof(double));
    read_csv_triplets("./test/A.csv",aI,aJ,aV);   
    printf("Loaded matrix A of size %i,%i with nnz: %i \n",m,n,nnzA);

    //read rhs.csv 
    int nnz = 0;
    nnz = read_size("./test/rhs.csv");
    double* rhs = calloc(nnz,sizeof(double));
    read_vector("./test/rhs.csv",rhs);
    printf("Loaded vector rhs of size %i \n",nnz);
    
    //Now call the solver
    double mu = 1;
    double delta = 1.e-10;
    double gamma = 1.e-10;
    
    double* sol = calloc(n+m,sizeof(double));
    build_and_solve_linear_system(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, sol, rhs);

    //Save the solution 
    write_vector_to_csv("./test/calculated_sol.csv",sol,n+m);
    printf("Saving calculated solution of size %i to csv\n",n+m); 
    //read solution.csv
    nnz = 0;
    nnz = read_size("./test/solution.csv");
    double* sol_mat = calloc(nnz,sizeof(double));
    read_vector("./test/rhs.csv",sol_mat); 
    printf("Loaded stored solution of size %i \n",nnz);
    
    double naccum = 0;
    int i = 0;
    double t = 0;
    for(i=0;i<nnz;i++)
    {
        t = (sol_mat[i]-sol[i]);
        naccum += (t*t);
    }
    naccum = sqrt(naccum);
    printf("Difference between saved and calculated sol %lf\n",naccum);    

}

int main (void)
{
    test_load_csv_and_solve_system();
}
