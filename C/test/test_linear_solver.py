#Try to call the library with ctypes
from ctypes import *
library = cdll.LoadLibrary('../Lib/liblinear_solvers.dylib');

nnz = 12
n   = 5

cnnz = c_int(nnz);
cn   = c_int(n);

nnzInt = c_int*nnz;
nnzDou = c_double*nnz;
nDou   = c_double*n;

x = create_string_buffer(n*8);

b  = nDou(8, 45, -3, 3, 19)
I = nnzInt(0,1,0,2,4,1,2,3,4,2,1,4)
J = nnzInt(0,0,1,1,1,2,2,2,2,3,4,4)
V = nnzDou(2.000000,3.000000,3.000000,-1.000000,4.000000,4.000000,-3.000000,1.000000,2.000000,2.000000,6.000000,1.000000)

library.solve_linear_system(x,I,J,V,cnnz,b,cn);
x_o = nDou.from_buffer(x);
print [i for i in x_o]

