#include "linear_solvers.h"
#include "umfpack.h"
#include <stdio.h>
#include "nscs_sp.h"
#include <OpenBlas/cblas.h>
#include "smatvec.h"
#include "test_util.h"


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
     int* Ai=(int*)calloc(nnz,sizeof(int));
     int* Ap=(int*)calloc(n+1,sizeof(int));
     double* Ax=(double*)calloc(nnz,sizeof(double)); 
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
void form_kkt_system_triplets( double mu,\
                              int* hI,\
                              int* hJ,\
                              double* hV,\
                              int nnzH,\
                              int* aI,\
                              int* aJ,\
                              double* aV,\
                              int nnzA,\
                              int m,\
                              int n,\
                              double delta,\
                              double gamma,\
                              int* I,\
                              int* J,\
                              double* V)
{ 

    int nnz;
    nnz     = 0; //Use nnz as a counter of the copied values 
    int i   = 0;

    //Set delta I
    for(i=0;i<m;i++)
    {
        I[nnz] = i;
        J[nnz] = i;
        V[nnz] = delta;
        nnz++;
    }
   
    //Copy  A'
    for(i=0;i<nnzA;i++)
    {
        I[nnz] = aJ[i]+m;
        J[nnz] = aI[i];
        V[nnz] = aV[i];
        nnz++;
    }
    
    //Copy A
    for(i=0;i<nnzA;i++)
    {
        I[nnz] = aI[i];
        J[nnz] = aJ[i]+m;
        V[nnz] = aV[i];
        nnz++;
    }
    
    //Copy -muH-gI
    //Since H is spd there are non zeros in the diagonal allready
    for(i=0;i<nnzH;i++)
    {
        I[nnz] = hJ[i]+m;
        J[nnz] = hI[i]+m;
        V[nnz] = -hV[i]*mu + ((hJ[i] == hI[i])?-gamma:0); //And add the diagonal
        nnz++;
    }

}

//This function forms the matrix to factorize and stores it in cco form in Ai, Aj, Av
//The vectors Ai,Aj,Av are allocated inside
int form_kkt_system( double mu,\
                     int* hI,\
                     int* hJ,\
                     double* hV,\
                     int nnzH,\
                     int* aI,\
                     int* aJ,\
                     double* aV,\
                     int nnzA,\
                     int m,\
                     int n,\
                     double delta,\
                     double gamma,\
                     int** Ai,\
                     int** Ap,\
                     double** Av)
{
    //[delta I          A ]
    //[A'      -muH-gammaI]

    //Allocate space for the new triplet form
    int nnz = m + 2*nnzA + nnzH;
    int* I = (int*)calloc(nnz,sizeof(int));
    int* J = (int*)calloc(nnz,sizeof(int));
    double* V = (double*)calloc(nnz,sizeof(double));
    
     form_kkt_system_triplets(mu, hI, hJ, hV, nnzH, aI, aJ, aV, nnzA, m, n, delta, gamma, I, J, V);
         
     //Convert to compressed column form
     //Allocate the space 
     (*Ai)=(csi*)calloc(nnz,sizeof(csi));
     (*Ap)=(csi*)calloc(n+m+1,sizeof(csi));
     (*Av)=(double*)calloc(nnz,sizeof(double)); 
     int* Map=NULL;

    //Execute the call to the compression routine
     int status  =  umfpack_di_triplet_to_col(n+m,n+m,nnz,I,J,V,*Ap,*Ai,*Av,Map);
     if(status!=0)
     {
        fprintf(stderr,"Pack COO to CSR failed %i\n",status);
     }
     //Clear the temporary triplet form 
     free(I);
     free(J);
     free(V);
     
     return status;
      
}

//This function receives the matrix to factor in cco form and retunrs
//the numeric factorization in Numeric
int factor_kkt_system(void** Numeric, int* Ai, int* Ap, double* Av, int n, int m)
{
 
     int status;
     //Now call umfpack to solve
     double *null = (double *) NULL ;
     //Define some stuff 
     void *Symbolic;
     //Do the symbolic analysis
     status = umfpack_di_symbolic (n+m, n+m, Ap, Ai, Av, &Symbolic, null, null) ;
     
     if(status != 0)
     {
        fprintf(stderr,"umfpack symbolic failed %i\n",status);
        if(Symbolic != NULL)
            umfpack_di_free_symbolic (&Symbolic) ;

         free(Ai);
         free(Ap);
         free(Av);
        return status;
     }
     else //symbolic succeeded
     { 
         //Do the numeric analysis
         status = umfpack_di_numeric (Ap, Ai, Av, Symbolic, Numeric, NULL,NULL);
         if(status != 0) //If there was an error free the stuff
         { 

            fprintf(stderr,"umfpack numeric failed %i\n",status); 
             umfpack_di_free_symbolic (&Symbolic) ;
             if(Numeric != NULL)
                umfpack_di_free_numeric (Numeric) ;
             free(Ai);
             free(Ap);
             free(Av);
             return status;
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
    status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Av, x, rhs, Numeric, NULL, NULL) ; 
    if(status != 0)
        fprintf(stderr,"umfpack solve failed %i\n",status);
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
 int loop(double *x, double* b, int n)
 {
     int i = 0;
     for(i=0;i<n;i++)
         x[i] = b[i];
 }


//Forms and solves the KKT system 
/**
 * This function forms the regularized KKT system and uses UMFPACK to solve it.
 *   [dI   A -b      ] dy   r1
 *   [-A' gI  c -I   ] dx   r2
 *   [b' -c'     -1 ] dt = r3
 *   [    mH   I    ] ds   r4
 *   [       k    t ] dk   r5
 *
 * @param mu
 * @param H the hessian of the barrier at the present primal point
 * @param A the constraint matrix
 * @param b the right hand side of the equality constraint
 * @param c the objective function 
 * @param tau (t) in the matrix above 
 * @param kappa (k) in the matrix above
 * @param delta (d) above
 * @param gamma (g) above
 * @param r1
 * @param r2
 * @param r3
 * @param r4
 * @param r5
 * @param dy pointer to a vec structure with pre-allocated space for m variables
 * @param dx 
 * @param dt 
 * @param ds 
 * @param dk 
 */
int solve_kkt_system(double mu,\
                     spmat H,\
                     spmat A,\
                     vec b,\
                     vec c,\
                     double tau,\
                     double kappa,\
                     double delta,\
                     double gamma,\
                     vec r1,\
                     vec r2,\
                     double r3,\
                     vec r4,\
                     double r5,\
                     vec dy,\
                     vec dx,\
                     double* dt,\
                     vec ds,\
                     double* dk )
{ 
   // Solves  [    A -b      ] dy   r1
   //         [-A'    c -I   ] dx   r2
   //         [b' -c'     -1 ] dt = r3
   //         [    mH   I    ] ds   r4
   //         [       k    t ] dk   r5


    //
    //Build the KKT system    
    //
    //    %Forms 
   //         [    A -b      ] dy   r1
   //         [-A'    c -I   ] dx   r2
   //         [b' -c'     -1 ] dt = r3
   //         [    mH   I    ] ds   r4
   //         [       h   1  ] dk   r6
   // with h = k/t r6 = r5/tau

   double h  = kappa/tau;
   //Scale by 1/kappa
   double r6 = r5/tau; 

    //Forms 
    //        [    A  -b      ] dy     r1
    //        [A' -mH -c      ] dx     r7
    //        [b' -c'  h      ] dt =   r8
    //        [   mH     I    ] ds     r4
    //        [        h   1  ] dk     r6
    // with r7 = -r2-r4; r8 = r3+r6;
    
    
    //void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
    
    //Allocate a working vector
    vec r7;
    r7.n = r2.n;
    r7.v = (double*)calloc(r7.n,sizeof(double));


    cblas_dcopy(r2.n,r2.v,1,r7.v,1);
    //Add r2 and r4 in r7
    cblas_daxpy(r2.n,1.,r4.v,1,r7.v,1);
    cblas_dscal(r7.n,-1,r7.v,1);  //#TODO: free r2
    //r7 = -(r2+r4);
    
    //r8 = r3+r6;
    double r8 = r3+r6; 
    //  
    //        [    A  -b     ] dy     r1
    //        [A' -mH -c     ] dx     r7
    //        [       h_2    ] dt =   r9
    //        [   mH     I   ] ds     r4
    //        [       h    1 ] dk     r6
    //  
    //[b' -c'] [     A]^{-1}[-b]     =  h_1b
    //         [A' -mH]     [-c] 

    //[b' -c'] [     A]^{-1}[r1]     =  r_7b
    //         [A' -mH]     [r7] 

    //We solve 
    // tm_1 =  [     A]^{-1}[b ]    
    //         [A' -mH]     [-c] 
    // and keep it to calculate h_1b and r_7b more efficiently
 
    //Call the solver 
    
    //Build the rhs
    double* rhs  = (double*)calloc(A.n+A.m,sizeof(double));
    int i = 0;
    for(i=0;i<A.m;i++)
        rhs[i] = b.v[i];
    for(i=0;i<A.n;i++)
        rhs[i+A.m] = -c.v[i];

    //Allocate space for the coo matrix
    csi* Ki=NULL;
    csi* Kp=NULL;
    double* Kv=NULL;
    //Calculate the final number of non zeros
    csi nnzK = A.m + 2*A.nnz + H.nnz;
    //Pointer for the factorization
    void* Numeric;
    //Pointer to the solution
    double* tm_1 = (double*)calloc(A.n+A.m,sizeof(double));
   
    int res;
    //From the kkt system
    res = form_kkt_system(mu, H.I, H.J, H.V, H.nnz, A.I, A.J, A.V, A.nnz, A.m, A.n, delta, gamma, &Ki, &Kp, &Kv);
    if(res!=0)
    {
        fprintf(stderr,"Error forming kkt system, %i \n",res);
        return INTERNAL_ERROR;
    }

    //Factor 
    res = factor_kkt_system(&Numeric, Ki, Kp, Kv, A.n, A.m);
    if(res!=0)
    {
        fprintf(stderr,"Error factoring kkt system %i, nnzk %i\n",res,Kp[nnzK]);
        return INTERNAL_ERROR;
    }
   
    //Solve
    res = solve_factored_system(Numeric, Ki, Kp, Kv, rhs, tm_1);
    if(res!=0)
    {
        fprintf(stderr,"Error solving kkt system %i \n",res);
        return INTERNAL_ERROR;
    }


    //Now calculate 
    // h_1b   = tm_1'*[-b;-c];
    double h_1b =  cblas_ddot(A.m,tm_1    ,1  ,b.v,1); 
    h_1b        += cblas_ddot(A.n,tm_1+A.m,1  ,c.v,1);
    h_1b        = -h_1b;
    
    // r_7b   = tm_1'*[r1;r7];
    double r_7b =  cblas_ddot(A.m,tm_1,1,r1.v,1); 
    r_7b        += cblas_ddot(A.n,tm_1+A.m,1,r7.v,1);

    free(r7.v);
 
    // r9  = r8 - r_7b;
    // h_2 = h-h_1b; 
    double r9   = r8-r_7b;
    double h_2  = h-h_1b;
   
    // %We now start the back substitution
    // dt  = r9/h_2;

    *dt  = r9/h_2;
   
    // %reuse tm_1
    // tm_1 = slv_aug([r1;r7]-dt*[-b;-c]);
    
    //Construct the rhs
    //put [r1;r7] into rhs
    cblas_dcopy(A.m,r1.v,1,rhs,1);
    cblas_dcopy(A.n,r7.v,1,rhs+A.m,1);

    //put r1+dt*b into rhs, and then r7+dt*c
    cblas_daxpy(A.m,*dt,b.v,1,rhs,1);
    cblas_daxpy(A.n,*dt,c.v,1,rhs+A.m,1);

    //Call the solver
    res = solve_factored_system(Numeric,Ki,Kp,Kv,rhs,tm_1);
    if(res!=0)
    {
        fprintf(stderr,"Error solving kkt system %i \n",res);
        return INTERNAL_ERROR;
    }

    //We dont need rhs anymore
    free(rhs);

    //We can fee Ki,Kp,Kv  and Numeric
    free(Ki);
    free(Kp);
    free(Kv); 
    umfpack_di_free_numeric(&Numeric);

    // %Extract dy dx from dt
    // dy   = tm_1(1:m);
    // dx   = tm_1(m+1:m+n);

    cblas_dcopy(A.m,tm_1,1,dy.v,1);
    cblas_dcopy(A.n,tm_1+A.m,1,dx.v,1);

    //We need the csr form of H for the product
    csi* Hi =(csi*)calloc(H.nnz,sizeof(csi));
    csi* Hp =(csi*)calloc(H.n+1,sizeof(csi)); //XXX: With free variables this can be smaller
    double* Hv = (double*)calloc(H.nnz,sizeof(double));  
    res = umfpack_di_triplet_to_col(H.m,H.n,H.nnz,H.I,H.J,H.V,Hp,Hi,Hv,NULL);     
    if(res!=0)
    {
        fprintf(stderr,"Error converting triplet form H to CRS %i \n",res);
        return INTERNAL_ERROR;
    }


    // %Now back substitute to get ds and dkappa
    // ds   = r4-mu*H*dx;
    // dk   = r6-h*dt;
    
    //Use tm_1(1:n) //XXX: With free variables this is no longer n
    //Calculate -mu*H*dx
    //clear tm_1
    for(i=0;i<H.m;i++)
        tm_1[i] = 0;

    dspmv(H.m,H.n,-mu,Hp,Hi,Hv,tm_1,dx.v); 
    //void dspmv(int m, int n, double alpha, int* pA, int * iA, double* vA, double* y, double* x)
    cblas_daxpy(H.n,1.,r4.v,1,tm_1,1);
    //Now ds is in tm_1;
    cblas_dcopy(H.m,tm_1,1,ds.v,1);
    *dk     = r6-h*(*dt);
    
    //Free the csr form of H
    free(Hi);
    free(Hp);
    free(Hv);
    //Free the work vector
    free(tm_1);

    // %Check the residuals
    // n_res_1 = norm(A*dx-dt*b-r1);
    // n_res_2 = norm(-A'*dy + dt*c -ds -r2);
    // n_res_3 = norm(b'*dy-c'*dx -dk -r3); 
    // n_res_4 = norm(mu*H*dx+ds-r4);
    // n_res_5 = norm(kappa*dt+tau*dk-r5);

  
    return 0;
}

//Copies In[0] to In[n-1] to Out[1] to Out[n-1]
//This is usefull to debug MATLAB's bugs in callib
void dummy_copy(double* Out, double* In, int n)
{
    cblas_dcopy(n,In,1,Out,1);
}
