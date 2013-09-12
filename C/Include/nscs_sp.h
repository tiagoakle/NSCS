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

/**
 * Structure to hold a dense vector
 */
typedef struct
{
    //Size 
    csi n;
    //Data
    double *v;
} vec;

//Allocates a vector of size n
vec* calloc_vec(int n);

//Allocates a matrix of size nnz
spmat* calloc_spmat(int m, int n, int nnz);
//Frees the vector 
void free_vec(vec v);
//Frees the arrays that make the matrix
void free_spmat(spmat A);

#endif
