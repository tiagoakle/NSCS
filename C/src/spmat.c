#include "spmat.h"
/**
 * Contains the methods used to create delete and operate
 * on spmatrix structures and vec structures
 */

//Allocates a matrix of size nnz
spmat* calloc_spmat(int m, int n, int nnz)
{
    spmat* mat = (spmat*)calloc(1,sizeof(spmat));
    if(mat == NULL) return NULL;

    mat->m = m;
    mat->n = n;
    mat->nnz = nnz;
    return mat;
}

/**
 * Frees a matrix stored in a spmat structure
 */
void free_spmat(spmat A)
{
    free(A.I);
    free(A.J);
    free(A.V);
}

