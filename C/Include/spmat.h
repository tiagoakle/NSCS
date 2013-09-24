#ifndef H_NSCS_SP
#define H_NSCS_SP

#include <common.h>
/*
 * This file contains the definitions for the
 * sructures that we use to hold sparse marices and vectors
 */

/**
 * Definition of a sparse matrix
 */
typedef struct 
{
   //Size of the matrix 
   csi m;
   csi n;
   csi nnz;
   //Matrix data
   csi* I;
   csi* J;
   double *V;

} spmat;


//Allocates a matrix of size nnz
int calloc_spmat(spmat* spmat, csi m, csi n, csi nnz);
//Frees the arrays that make the matrix
void free_spmat(spmat A);

#endif
