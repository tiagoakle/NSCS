#include "linear_solvers.h"
#include "umfpack.h"
#include <stdio.h>

//This function builds the compressed column row structure and 
//    then calls umfpack to solve. Ax=b
//    The parameters are:
//        x a pointer for the resulting vector
//        I,J,V vectors with coo format of A
//        b     vector of values of b
//        n     size of the system
//
int solve_linear_system(double* x, int* I, int* J, double* V, int nnz, double *b ,int n)
{
  
     //Convert to compressed column form
     //Allocate the space 
     int* Ai=calloc(nnz,sizeof(int));
     int* Ap=calloc(n+1,sizeof(int));
     double* Ax=calloc(nnz,sizeof(double)); 
     int* Map=NULL;
     //Execute the call
     int status;

     status  =  umfpack_di_triplet_to_col(n,n,nnz,I,J,V,Ap,Ai,Ax,Map);
     if(status != 0)
        return status;
     
     //Now call umfpack to solve
     double *null = (double *) NULL ;
     //Define some stuff 
     void *Symbolic, *Numeric ;
     //Do the symbolic analysis
     status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
     if(status != 0)
     {
        if(Symbolic != NULL)
            umfpack_di_free_symbolic (&Symbolic) ;

         free(Ai);
         free(Ap);
         free(Ax);
        return status;
     }
     else //symbolic succeeded
     { 
         //Do the numeric analysis
         status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL,NULL) ;
         if(status != 0) //If there was an error free the stuff
         { 
             umfpack_di_free_symbolic (&Symbolic) ;
             if(Numeric != NULL)
                umfpack_di_free_numeric (&Numeric) ;
             free(Ai);
             free(Ap);
             free(Ax);
             return status;
         }
         else //Numeric succeeded 
         { 
            //Delete the structure that holds the symbolic analysis
            umfpack_di_free_symbolic (&Symbolic) ;
            // Now solve the system
            status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL) ;
            umfpack_di_free_numeric (&Numeric) ;
         }    
     }
     //Free the compressed column form 
     free(Ai);
     free(Ap);
     free(Ax);
     return status;
    
}

//Forms the KKT system in the allocated working vectors I,J,V 
int form_kkt_system_triplets( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int* I, int* J, double* V)
{ 

    int nnz;
    nnz     = 0; //Use nnz as a counter of the copied values 
    int i   = 0;

    //Allocate delta I
    for(i=0;i<m;i++)
    {
        I[nnz] = i;
        J[nnz] = i;
        V[nnz] = delta;
        nnz++;
    }
   
    //Allocate the A'
    for(i=0;i<nnzA;i++)
    {
        I[nnz] = aJ[i]+m;
        J[nnz] = aI[i];
        V[nnz] = aV[i];
        nnz++;
    }
    
    //Allocate the A
    for(i=0;i<nnzA;i++)
    {
        I[nnz] = aI[i];
        J[nnz] = aJ[i]+m;
        V[nnz] = aV[i];
        nnz++;
    }
    
    //Allocate the -muH 
    //Since H is spd there are non zeros in the diagonal allready
    for(i=0;i<nnzH;i++)
    {
        I[nnz] = hJ[i]+m;
        J[nnz] = hI[i]+m;
        V[nnz] = -hV[i]*mu + ((hJ[i] == hI[i])?-gamma:0); //And add the diagonal
        nnz++;
    }

    return 0;

}

//This function forms the matrix to factorize and stores it in cco form in Ai, Aj, Av
int form_kkt_system( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, int** Ai, int** Ap, double** Av)
{
    //[delta I          A ]
    //[A'      -muH-gammaI]

    //Allocate space for the new triplet form
    int nnz = m + 2*nnzA + nnzH;
    int* I = calloc(nnz,sizeof(int));
    int* J = calloc(nnz,sizeof(int));
    double* V = calloc(nnz,sizeof(double));
    
    int status = form_kkt_system_triplets(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, I, J, V);
     
     //Convert to compressed column form
     //Allocate the space 
     (*Ai)=calloc(nnz,sizeof(int));
     (*Ap)=calloc(n+m+1,sizeof(int));
     (*Av)=calloc(nnz,sizeof(double)); 
     int* Map=NULL;

    //Execute the call to the compression routine
     status  =  umfpack_di_triplet_to_col(n+m,n+m,nnz,I,J,V,*Ap,*Ai,*Av,Map);

     //Clear the temporary triplet form 
     free(I);
     free(J);
     free(V);
     
     return status;
      
}

//This function receives the matrix to factor in cco form and retunrs
//the numeric factorization in Numeric
int factor_kkt_system(void** Numeric, int* Ai, int* Ap, double* Av, int n)
{

 
     int status;
     //Now call umfpack to solve
     double *null = (double *) NULL ;
     //Define some stuff 
     void *Symbolic;
     //Do the symbolic analysis
     status = umfpack_di_symbolic (n, n, Ap, Ai, Av, &Symbolic, null, null) ;

     printf("Symbolic  KKT SYSTEM %i\n",status);     
     if(status != 0)
     {
        if(Symbolic != NULL)
            umfpack_di_free_symbolic (&Symbolic) ;

         free(Ai);
         free(Ap);
         free(Av);
        return status - 100;
     }
     else //symbolic succeeded
     { 
         //Do the numeric analysis
         status = umfpack_di_numeric (Ap, Ai, Av, Symbolic, Numeric, NULL,NULL) ;

         printf("Numeric  KKT SYSTEM %i\n",status);     
         if(status != 0) //If there was an error free the stuff
         { 
             umfpack_di_free_symbolic (&Symbolic) ;
             if(Numeric != NULL)
                umfpack_di_free_numeric (Numeric) ;
             free(Ai);
             free(Ap);
             free(Av);
             return status - 200;
         }
         else //Numeric succeeded 
         { 
            //Delete the structure that holds the symbolic analysis
            umfpack_di_free_symbolic (&Symbolic) ;
         }    
     }
     return status;  
}

//Solves the linear system using the factorization
int solve_factored_system(void* Numeric, int* Ai, int* Ap, double* Av, double* rhs, double* x)
{
    int status; 
    // Now solve the system
    status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Av, x, rhs, Numeric, NULL, NULL) ; 
   
    return status; 
}

//Clear the memory used for the factored matrix
int free_factorization(void* Numeric, int* Ai, int* Ap, double* Av)
{
     
     umfpack_di_free_numeric (&Numeric) ;
     //Free the compressed column form 
     free(Ai);
     free(Ap);
     free(Av);
 
}

//This is for debugging
int build_and_solve_linear_system( double mu, int* hI, int* hJ, double* hV, int nnzH, int* aI, int* aJ, double* aV, int nnzA, int m, int n, double delta, double gamma, double* sol, double * rhs)
{
    int* Ai = NULL;
    int* Ap = NULL;
    double* Av = NULL;
    void* Numeric = NULL;

     int ret = form_kkt_system(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, &Ai, &Ap, &Av);
  //   int nnz = m + 2*nnzA+nnzH;
  //       ret = solve_linear_system(x, Ai, Aj, Av, nnz, rhs ,n +m)
    printf("Form KKT SYSTEM %i\n",ret);     
     if (ret != 0)
        sol[0] = -1000+ret;     
     ret = factor_kkt_system(&Numeric, Ai, Ap, Av, n+m);
    printf("Factor KKT SYSTEM %i\n",ret);     
     if (ret != 0)
        sol[1] = -2000+ret;
     ret = solve_factored_system(Numeric, Ai, Ap, Av, rhs, sol);
    printf("Solve KKT SYSTEM %i\n",ret);     
     if (ret != 0)
        sol[2] = -3000+ret;
     //ret = free_factorization(Numeric, Ai,Ap,Av);
    return ret;
}

//This is for debugging
int loop(double *x, double* b, int n)
{
    int i = 0;
    for(i=0;i<n;i++)
        x[i] = b[i];
}
