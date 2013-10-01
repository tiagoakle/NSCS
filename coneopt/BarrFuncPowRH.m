function [F,G,H,feas] = BarrFuncPowRH(x1,x2,x3,K,compute)

% the function 
%  F = Einsein-Hessian barrier from Roland Hildebrand
%  for the product of k cones abs(x1) <= x2^alpha x3^(1-alpha) 
% (each of x1, x2, x3 and alpha is a row vector of size k)
%
%  output:
%   F    = function value
%   G    = Elements of gradient
%   H    = Elements of Hessian 
%   feas = 0 if point is feasible, -1 if not
%
%   Notice below the way that G and H are computed and output.
%   For example, H is output as a vector. This is because, it must later
%   be gathered with Hessians from other cones computed by other barrierfunctions
%   This we do using one big assembly-procedure where indices are pre-computed 
%   to make it reasonably fast. Therefore, G and H should just be output in 
%   the way seen below.  
%
%  input:
%   x1,x2,x3 are vectors of length k
%   so that ( x1(j),x2(j),x3(j) ) is in powercone( alpha(j) )
%   alpha is therefore also a vector of length k specifying the parameters
%   of each of the k powercones. We assume that alpha is passed through 
%   the struct K.
%
%   compute = vector with three elements
%    if compute(1) > 0, F is computed
%    if compute(2) > 0, G is computed
%    if compute(3) > 0, H is computed
%   this appears so that we can specify not to compute F/G/H when they are
%   not needed, e.g. if we only need to check feasibility in line search.
% 
%  First some things that can be computed outside this function
%  (I include it here so that the function works by itself. These
%  structures are normally passed into the function through the struct K). 

alpha = K.powc.alph;

% [F,F1,F2,s] = EHbarrier(x,y,z,p)
% computess barrier F, gradient/hessian F1, F2 and boolean s for
% feasibility for the cone |z| < (sqrt(p)*x)^(1/p)*(sqrt(q)*y)^(1/q)
% Need to equate to abs(x1) <= x2^alpha x3^(1-alpha) 
%                   abs(x1) <= (sqrt(p)*x2/sqrt(p))^alpha sqrt(q)*x3/sqrt(q))^(1-alpha) 
% hence choose x=x2*sqrt(alpha) ; y=x3*sqrt(1-alpha) ; z=x1 ; p=1/alpha

k = length(alpha);

feas = zeros(1, k);
F = zeros(1, k);
G = zeros(3, k);
H = zeros(9, k);
for i=1:k % slow cone by cone computation ...
   [f,F1,F2,s] = EHbarrierRH(x2(i)*sqrt(alpha(i)),x3(i)*sqrt(1-alpha(i)),x1(i),1/alpha(i));
   feas(1,i) = s;
   F(1,i) = f;
   % define scaling such that bar(v) = RH(scaling*x)
   scaling = [0 sqrt(alpha(i)) 0 ; 0 0 sqrt(1-alpha(i)) ; 1 0 0];
   G(:,i) = scaling'*F1;
   hess = scaling'*F2*scaling;
   H(:,i) = hess(:);
%   if s>0  % two simple check for grad and hessian  
%      hess*[x1(i);x2(i);x3(i)]+G(:,i) 
%      -[x1(i) x2(i) x3(i)]*G(:,i)
%   end
end
F = sum(F);
if all(feas)
    feas = 0;
else 
    feas = -1;
end
G=G(:);
H=H(:);
