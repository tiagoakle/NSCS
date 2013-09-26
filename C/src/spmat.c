#include "spmat.h"
#include "common.h"
#include <stdio.h>

/**
 * Contains the methods used to create and delete
 * spmat structures
 */

/**
 * Allocates a matrix of size m,n with nnz non zeros
 */
int calloc_spmat(spmat* spmat, csi m, csi n, csi nnz)
{
    if(spmat == NULL){fprintf(stderr,"Missing argument\n"); return MISSING_ARGUMENT;}
    spmat->m = m;
    spmat->n = n;
    spmat->nnz = nnz;
    spmat->I = (csi*)calloc(nnz,sizeof(csi));
    if(spmat->I==NULL){fprintf(stderr,"Out of memory\n");return OUT_OF_MEMORY;}
    spmat->J = (csi*)calloc(nnz,sizeof(csi));
    if(spmat->J==NULL){fprintf(stderr,"Out of memory\n");return OUT_OF_MEMORY;}
    spmat->V = (double*)calloc(nnz,sizeof(double));
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

