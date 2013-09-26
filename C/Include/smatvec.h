#include "common.h"
void dspmv(int n, int m, double alpha, int* pA, int * iA, double* vA, double* y, double *x);
void dsTpmv(int n, int m, double alpha, int* pA, int * iA, double* vA, double* y, double *x);
void dspmvcoo(csi nnz, double alpha, int* Ai, int * Aj, double* Av, double* y, double* x);
void dsptmvcoo(csi nnz, double alpha, int* Ai, int * Aj, double* Av, double* y, double* x);
