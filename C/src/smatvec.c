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
