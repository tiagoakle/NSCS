#include "smatvec.h"
#include "common.h"
#include "stdio.h"
/**
 * Naive sparse matrix dense vector multiply
 * y = alpha Ax + y
 */
void dspmv(int m, int n, double alpha, int* pA, int * iA, double* vA, double* y, double* x)
{
    //Direct product
    
    csi i = 0;
    csi c = 0;
    double xc = 0;
       
    //Iterate over the columns
    for(c=0;c<n;c++)
    {
        xc = x[c];
        while(i<pA[c+1]) //Jump zero columns
        {
            y[iA[i]] += alpha*vA[i]*xc;
            i++;
        }
    }

}

/**
 *  Naive implementation
 *  of the transpose product.
 *
 *  y = alpha A'x + y
 *
 */
void dsTpmv(int m, int n, double alpha, int* pA, int * iA, double* vA, double* y, double* x)
{
    //Transpose product
    csi i = 0;
    csi c = 0;
    
    //Iterate over the columns
    for(c=0;c<n;c++)
    {
        while(i<pA[c+1]) //Jump zero columns
        {
            y[c] += alpha*vA[i]*x[iA[i]];
            i++;
        }
    }
    
}

/**
 * Naive sparse matrix dense vector multiply in coordinate form
 * y = alpha Ax + y
 */
void dspmvcoo(csi nnz, double alpha, int* Ai, int * Aj, double* Av, double* y, double* x)
{
    //Direct product 
    csi i = 0; 
    //Iterate over all non zeros
    for(i=0;i<nnz;i++)
    {
            y[Ai[i]] += alpha*Av[i]*x[Aj[i]];
    }

}

/**
 * Naive sparse Transpose matrix dense vector multiply in coordinate form
 * y = alpha A'x + y
 */
void dsptmvcoo(csi nnz, double alpha, int* Ai, int * Aj, double* Av, double* y, double* x)
{
    //Direct product
    
    csi i = 0;
       
    //Iterate over all non zeros
    for(i=0;i<nnz;i++)
    {
            y[Aj[i]] += alpha*Av[i]*x[Ai[i]];
    }

}


