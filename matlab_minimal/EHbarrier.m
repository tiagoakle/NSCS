function [F,F1,F2,s] = EHbarrier(x,y,z,p)
% function [F,F1,F2,s] = EHbarrier(x,y,z,p)
%
% computes the Einstein-Hessian barrier on the power cone with parameter p
% the cone is given by
% |z| < (sqrt(p)*x)^(1/p)*(sqrt(q)*y)^(1/q)
% s is true if the arguments lie in the interior of the cone
% F is the barrier value
% F1 is the gradient as a column vector
% F2 is the Hessian
% the barrier has parameter 3s/(s+1), where s = max(p,q)
q = 1/(1-1/p);
if p < 2,
    r = p;
    p = q;
    q = r;
end,
F = 0;
F1 = zeros(3,1);
F2 = zeros(3);
if min([x,y]) <= 0,
    s = false;
    return;
end,
t = (sqrt(p)*x)^(-1/p)*(sqrt(q)*y)^(-1/q)*z;
if abs(t) >= 1,
    s = false;
    return;
end,
s = true;
if t == 0,
    phi = 1/2*((p+1)/p*log(p+1)+(q+1)/q*log(q+1));
    F = (-(p+1)/p*log(sqrt(p)*x) - (q+1)/q*log(sqrt(q)*y) + phi)*p/(p+1);
    F1 = [-(p+1)/(p*x); -(q+1)/(q*y); 0]*p/(p+1);
    F2 = [(p+1)/(p*x^2), 0, 0; 0, (q+1)/(q*y^2), 0; 0, 0, exp(2*phi)/(2*p+2*q+1)/(sqrt(p)*x)^(2/p)/(sqrt(q)*y)^(2/q)]*p/(p+1);
else
    if abs(t) < 10^(-4),
        k1 = (p+1)^(1/p)*(q+1)^(1/q)*t^2;
        k = k1 + k1^2*(1/(p*(p+1))+1/(q*(q+1)));
    else
        k1i = (p+1)^(-1/p)*(q+1)^(-1/q)*t^(-2);
        ki = max([k1i - 1/(p*(p+1))+1/(q*(q+1)), (k1i*(p+1)^(1/p)*(q+1)^(1/q)-1)/3]);
        f = (ki + 1/(p+1))^(1/p)*(ki + 1/(q+1))^(1/q);
        while abs(f-k1i)/k1i > 10^(-14),
            f1 = f/p/(ki + 1/(p+1)) + f/q/(ki + 1/(q+1));
            ki = ki + (k1i-f)/f1;
            f = (ki + 1/(p+1))^(1/p)*(ki + 1/(q+1))^(1/q);
        end,
        k = 1/ki;
    end,
    phi = 1/2*((p+1)/p*log(k+p+1)+(q+1)/q*log(k+q+1));
    t2p = 2*k*(k+p+1)*(k+q+1)/(3*k+2*(p+q)+1)-k;
    M = diag([1/(p*x), 1/(q*y), 1/z]);
    F = (-(p+1)/p*log(sqrt(p)*x) - (q+1)/q*log(sqrt(q)*y) + phi)*p/(p+1);
    F1 = M*[-(p+1+k); -(q+1+k); k]*p/(p+1);
    F2 = M*[(p+k)*(p+1)+t2p, k+t2p, -(k+t2p); k+t2p, (q+k)*(q+1)+t2p, -(k+t2p); -(k+t2p), -(k+t2p), t2p]*M*p/(p+1);
end,
