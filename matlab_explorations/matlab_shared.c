#include "matlab_shared.h"
#include <stdio.h>

int print_hello()
{
    printf("Hello\n");
}

//Load a triplet from the matlab memory
void TransposeInTripletForm(int* I, int* J, double* Vals, int nnz)
{
    if(nnz == 0)
        return;
    //Copy from J to I
    int t,i;
    for(i=0;i<nnz;i++)
    {
        t = I[i];
        I[i] = J[i];
        J[i] = t;
    }

}

//Copy the n values in the array that starts at in
//to the array that starts at out
void CopyVals(double *pOut, double* in, int n)
{
    int i = 0;
    for(i=0;i<n;i++)
        *(pOut++) = *(in++);
        
}

