#include <stdio.h>
#include <stdlib.h>
#include "SuiteSparse/CSparse/Include/cs.h"

//Generates a random matrix with nnz non zeros in 
//triplet form readable by suite sparse
cs *create_triplet_matrix(csi m, csi n, csi nzmax)
{
    //Allocate the matrix 
   csi values = 1; //Indicate we want to allocate space for the non zeros
   csi triplet= 1; //Indicate we want to store in triplet form
   cs *T = cs_spalloc(m,n,nzmax,values,triplet);
   int i = 0;
   for(i=0;i++;i<nzmax)
   {
        csi ix_i = (csi)((rand()/(double) (RAND_MAX+1))*n); 
        csi ix_j = (csi)((rand()/(double) (RAND_MAX+1))*m);
        double val = rand()/(double) (RAND_MAX+1);
        cs_entry(T,ix_i,ix_i,val);
   }
   return T;
}


int main(int argn, char* args)
{
   cs *T, *C;
   //Call the test functions 
   printf("Will crete the matrix\n");
   T = create_triplet_matrix(10,10,90);
   printf("Created the matrix\n");
   //Convert the matrix to compressed column form
   C = cs_compress(T);
   printf("Compressed the matrix\n");
   //Free the matrix
   T = cs_spfree(T);
   C = cs_spfree(C);
   printf("Freed the matrix\n");
}
