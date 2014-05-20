%Symbolic calculation of the Cholesky factor of the hessian
%We will have to calculate sqrt(vH^{-1}(x)v) to measure the centrality 
%so we should use the symbolic toolbox for the cholesky factorization of H

syms x y z s1 s2 s3 g1 g2 g3 mu real
f = -log(z*log(y/z)-x)-log(y)-log(z);
g = gradient(f)
s = [s1;s2;s3];
H = hessian(f);
Hi = inv(H);
v = (s/mu+g)'*Hi*(s/mu+g)
v = simplify(v,2000)
