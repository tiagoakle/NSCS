#include <stdio.h>
#include "linear_solvers.h"

//Test to call linear_solvers
int n = 5 ;
int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
int Aj [ ] = { 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4} ;
int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double b [ ] = {8., 45., -3., 3., 19.} ;
double x [5] ;

int main (void)
{
    solve_linear_system(x, Ai, Aj, Ax, 12, b , n);
    size_t i = 0;
    printf("RHS \n");
    for(i=0;i<n;i++)
        printf("%f\n",b[i]);

    printf("Solution \n");
    for(i=0;i<n;i++)
        printf("%f\n",x[i]);
}
