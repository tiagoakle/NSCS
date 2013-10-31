%Ecos test
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


