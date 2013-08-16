#include "linear_solvers.h"
#include "umfpack.h"
#include <stdio.h>

//This is a test function that mutates the n entries in x by copying
//from b, and mutates the nnz entries in I and V by adding one to them
void mutate_mem(double* x, int *I, double *V, double *b, int nnz, int n)
{
    int i = 0;
    for(i=0;i<n;i++)
        *(x++) = *(b++);
    for(i=0;i<n;i++)
    {
        (*(I++))++;
        (*(V++))++;
    }
}

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
  
 //Debug File 
    FILE *f;
    f = fopen("Data_dump.txt", "w");
    int i;
    fprintf(f,"\n---B -----");
    for(i=0;i<n;i++)
        fprintf(f, "%f,",b[i]);

    fprintf(f,"\n---I-----");
    for(i=0;i<nnz;i++)
        fprintf(f, "%i,",I[i]);

    fprintf(f,"\n---J-----");
    for(i=0;i<nnz;i++)
        fprintf(f, "%i,",J[i]);

    fprintf(f,"\n---V-----");
    for(i=0;i<nnz;i++)
        fprintf(f, "%f,",V[i]);
    fprintf(f,"\n nnz: %i \n n: %i ",nnz,n);
    fflush(f);

     //Convert to compressed column form
     //Allocate the space 
     int* Ai=calloc(nnz,sizeof(int));
     int* Ap=calloc(n+1,sizeof(int));
     double* Ax=calloc(nnz,sizeof(double)); 
     int* Map=NULL;
     //Execute the call
     int status;

     status  =  umfpack_di_triplet_to_col(n,n,nnz,I,J,V,Ap,Ai,Ax,Map);
     fprintf(f,"\nStatus to col %i",status);
     fflush(f);
    
    //Now call umfpack to solve
     double *null = (double *) NULL ;
     //Define some stuff 
     void *Symbolic, *Numeric ;
     //Do the symbolic analysis
     //
     //Desperate try copy all the data to new buffers
     int* Apnew = calloc(n+1,sizeof(int));
     int* Ainew = calloc(nnz,sizeof(int));
     double* Axnew = calloc(nnz,sizeof(double));
     double* bnew  = calloc(n, sizeof(double));
     double* xnew  = calloc(n,sizeof(double));

     for(i=0;i<n+1;i++)
        Apnew[i] = Ap[i];
     for(i=0;i<n;i++)
        bnew[i] = b[i];
     for(i=0;i<nnz;i++)
     {
        Ainew[i] = Ai[i];
        Axnew[i] = Ax[i];
     }

   
     //
     status = umfpack_di_symbolic (n, n, Apnew, Ainew, Axnew, &Symbolic, NULL, NULL) ;
     fprintf(f,"\nStatus sym %i",status);
     fflush(f);

    fprintf(f,"\n----Ai-----");
    for(i=0;i<nnz;i++)
        fprintf(f, "%i,",Ainew[i]);

    fprintf(f,"\n----Ax---");
    for(i=0;i<nnz;i++)
        fprintf(f, "%f,",Axnew[i]);

    fprintf(f,"\n----Ap----");
    for(i=0;i<n+1;i++)
        fprintf(f, "%i,",Apnew[i]);
    
    fflush(f); 


        //Do the numeric analysis
     status = umfpack_di_numeric (Apnew, Ainew, Axnew, Symbolic, &Numeric, NULL,NULL) ;
     fprintf(f,"\nStatus num %i",status);
     fflush(f);
   
     //Delete the structure that holds the symbolic analysis
     umfpack_di_free_symbolic (&Symbolic) ;
     // Now solve the system
     status = umfpack_di_solve (UMFPACK_A, Apnew, Ainew, Axnew, xnew, bnew, Numeric, null, null) ;
     
     fprintf(f,"\nStatus solve %i",status);
     fflush(f);
     
     umfpack_di_free_numeric (&Numeric) ;

     //Copy to x
    for(i=0;i<n;i++)
      x[i] = xnew[i];

    
       fprintf(f,"\n---X-----");
    for(i=0;i<n;i++)
        fprintf(f, "%f,",x[i]);

     free(Ai);
     free(Ap);
     free(Ax);
     free(Ainew);
     free(Axnew);
     free(Apnew);
     free(bnew);
     free(xnew);

   fclose(f);
}

