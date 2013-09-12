#include "nscs_sp.h"
/**
 * Contains the methods used to create delete and operate
 * on spmatrix structures and vec structures
 */

//Allocates a vector of size n
vec* calloc_vec(int n)
{
    vec* v = (vec*)calloc(1, sizeof(vec));
    if(v == NULL) return NULL;

    v->v   = (double*)calloc(n, sizeof(double));
    v->n   = n;
    return v;
}

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

void free_vec(vec v)
{
    free(v.v);
}
