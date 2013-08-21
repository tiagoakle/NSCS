[ret,warn] = loadlibrary('../lib/liblinear_solvers','../include/linear_solvers.h');

A       = sprand(n,n,1);
n       = 300;
tests   = 1000;
lpX     = libpointer('doublePtr',zeros(n,1));
rhs1    = ones(n/2,1);
rhs2    = zeros(n/2,1);
rhs     = [rhs1;rhs2];
ret     = calllib('liblinear_solvers','loop',lpX,rhs,n);
err     = norm(rhs-lpX.value)

rhs1    = ones(n/2,1);
rhs2    = zeros(n/2,1);
rhs     = [rhs1;rhs2];
ret     = calllib('liblinear_solvers','loop',lpX,rhs,n);
err     = norm(rhs-lpX.value)


unloadlibrary liblinear_solvers


