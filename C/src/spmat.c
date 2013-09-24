#include "spmat.h"
#include "common.h"
#include <stdio.h>

/**
 * Contains the methods used to create delete and operate
 * on spmatrix structures and vec structures
 */

//Allocates a matrix of size nnz
int calloc_spmat(spmat* spmat, csi m, csi n, csi nnz)
{
    if(spmat == NULL){fprintf(stderr,"Missing argument\n"); return MISSING_ARGUMENT;}
    spmat->m = m;
    spmat->n = n;
    spmat->nnz = nnz;
    spmat->I = calloc(nnz,sizeof(csi));
    if(spmat->I==NULL){fprintf(stderr,"Out of memory\n");return OUT_OF_MEMORY;}
    spmat->J = calloc(nnz,sizeof(csi));
    if(spmat->J==NULL){fprintf(stderr,"Out of memory\n");return OUT_OF_MEMORY;}
    spmat->V = calloc(nnz,sizeof(double));
    if(spmat->V==NULL){fprintf(stderr,"Out of memory\n");return OUT_OF_MEMORY;}
    return OK;
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

