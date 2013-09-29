#include <stdio.h>
#include "test_util.h"
#include "nscs.h"
#include "common.h"
#include "umfpack.h"

//Loads a problem , A,b,c from ./test_data/minentropy/...
//generated by build_minentropy_example_files and calls
//nscs to solve it
int main(void)
{
    int m,n,nnz;
    int* Ai;
    int* Aj;
    double* Av;
    double* b;
    double* c;
    double* x;
    double* s;
    double* y;
    double kappa;
    double tau;
    int* nK;
    int* tK;
    int k_count;
    bool wy,ws,wt,wk;
    double p_relstop, d_relstop, g_relstop;
    double relstops[3];
    //Read the dimensions of the problem
    read_matrix_bin_size("./test/test_data/minentropy/A.bin", &m, &n, &nnz);
    //Allocate the space
    Ai = (int*)calloc(nnz,sizeof(int));
    Aj = (int*)calloc(nnz,sizeof(int));
    Av = (double*)calloc(nnz,sizeof(double));
     c  = (double*)calloc(n,sizeof(double));
     b  = (double*)calloc(m,sizeof(double));
     x  = (double*)calloc(n,sizeof(double));
     s  = (double*)calloc(n,sizeof(double));
     y  = (double*)calloc(m,sizeof(double));
    
    //Sanity check on A 
     int* Map = NULL;
     int* P = (int*)calloc(n+1,sizeof(int));
     int* I = (int*)calloc(nnz,sizeof(int));
     double* V = (double*)calloc(nnz,sizeof(double));
     int status =  umfpack_di_triplet_to_col(m,n,nnz,Ai,Aj,Av,P,I,V,Map);
     printf("Status of compression %i\n",status);

    //Read the A,b,c,x
    read_matrix_bin("./test/test_data/minentropy/A.bin", Ai,Aj,Av);
    read_vector_bin("./test/test_data/minentropy/c.bin",c,n); 
    read_vector_bin("./test/test_data/minentropy/b.bin",b,m);
    read_vector_bin("./test/test_data/minentropy/x0.bin",x,n);
    //Read the relstop parameters
    read_vector_bin("./test/test_data/minentropy/relstop.bin",relstops,3);
    p_relstop = relstops[0];
    d_relstop = relstops[1];
    g_relstop = relstops[2];
   
    //Build the cone arrays
    tK = (int*)calloc(n/3,sizeof(int));
    nK = (int*)calloc(n/3,sizeof(int));
    int i;
    for(i=0;i<n/3;i++)
    {
        tK[i] = 3;
        nK[i] = 3;
    }

    problem_t prob;
    spmat A;
    A.I = Ai;
    A.J = Aj;
    A.V = Av;
    A.m = m;
    A.n = n;
    A.nnz = nnz;
    prob.A = A;
    prob.b = b;
    prob.c = c;
    prob.tK = tK;
    prob.nK = nK;
    prob.k_count = n/3;
    
    //Now build the paramters structure
    parameters_t pars;
    pars.print = true;
    pars.max_iter = 500;
    pars.max_center_iter = 50;
    pars.theta = 0.8;
    pars.lscaff  = 0.94;
    pars.lsccent = 0.5;
    pars.eta     = 0.9995;
    pars.max_backtrack = 300;
    pars.delta = 1.e-10;
    pars.beta  = 0.99;
    pars.beta      = 0.99;
    pars.p_relstop = p_relstop;
    pars.d_relstop = d_relstop;
    pars.g_relstop = g_relstop;
    pars.p_rho = 1e-6; 
    pars.d_rho = 1e-6;
    pars.a_rho = 1e-6; 
    pars.rho_g = 1e-6;
    pars.rhoI  = 1e-6;
    pars.rhoM  = 1e-6;
    
    nscs(&prob,&pars,y,x,&tau,s,&kappa,false,false,false,false);

}
