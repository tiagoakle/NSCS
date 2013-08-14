#include "SuiteSparse/CSparse/Include/cs.h"
#include "linear_solvers.h"
//This function builds the compressed column row structure and 
//    then calls umfpack to solve. Ax=b
//    The parameters are:
//        x a pointer for the resulting vector
//        I,J,V vectors with coo format of A
//        b     vector of values of b
//        n     size of the system
//
void solve_linear_system(double* x, int* I, int* J, double* V, double *b ,int n)
{
    int i =0;
}

