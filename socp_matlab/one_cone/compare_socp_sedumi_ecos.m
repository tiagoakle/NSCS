clear all
%Generates a problem with socp, then calls sedumi and ecos
one_cone_socp_test

%Now call sedumi 
K.q = n;
sedumi(A,b,c,K);

%Now call ecos
cvx_begin
    cvx_solver 'ecos'
    variable x(n)
    dual variable y
    dual variable s
    minimize(c'*x)
    subject to
    y: A*x == b
    s: x(1) >= norm(x(2:n))
  cvx_end

%Sanity check sedumi call
cvx_begin
    cvx_solver 'sedumi'
    variable x(n)
    dual variable y
    dual variable s
    minimize(c'*x)
    subject to
    y: A*x == b
    s: x(1) >= norm(x(2:n))
  cvx_end


G = -speye(n);
h = zeros(n,1);
dims.l = [];
dims.q = n;
ecos(c,G,h,dims,A,b);
