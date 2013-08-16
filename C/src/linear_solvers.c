#include "linear_solvers.h"
//This function builds the compressed column row structure and 
//    then calls umfpack to solve. Ax=b
//    The parameters are:
//        x a pointer for the resulting vector
//        I,J,V vectors with coo format of A
//        b     vector of values of b
//        n     size of the system
//
void solve_linear_system(double* x, int* I, int* J, double* V, int nnz, double *b ,int n)
{
    //For now just copy b to x
    size_t i;
    for(i=0;i<n;i++)
        *(x++) = *(b++);
    
    // cs* A = calloc(1,sizeof(cs)); 
    //XXX: This is doubling the memory ussage!
    //     we have I,J,V in the heap and all the space in cs_spalloc
    //Allocate a csparse structure
    cs *T = cs_spalloc(n,n,nnz,1,1);
    //Assign the entries
    i = 0;
    for(i=0;i<nnz;i++)
     cs_entry(T,I[i],J[i],V[i]);
    //Compress to column form
    cs *C = cs_compress(T);
    
    //Now call umfpack!
    double *null = (double *) NULL ;
    
    void *Symbolic, *Numeric ;
    (void) umfpack_dl_symbolic (n, n, C->p, C->i, C->x, &Symbolic, NULL, NULL) ;
    (void) umfpack_dl_numeric (C->p, C->i, C->x, Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    (void) umfpack_dl_solve (UMFPACK_A, C->p, C->i, C->x, x, b, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
}

