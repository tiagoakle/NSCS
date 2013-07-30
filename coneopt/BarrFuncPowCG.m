function [F,G,H,feas] = BarrFuncPowCG(x1,x2,x3,K,compute)

% the function 
%  F = -log( x2^(2*alpha)*x3^(2*(1-alpha)) - x1^2 )  -(1-alpha)*log(x2) -alpha*log(x3)
%  which has barrier param equal to 3 (Chares/Glineur).
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

k     = length(x1);
alpha = K.alpha;
tmp   = [alpha' ; 1-alpha'];
a1    = flipud(tmp);
a2    = 2*tmp;
ga1   = [ones(1,k);tmp]; %k is the number of power cones (i.e. the length of x1,x2,x3).
ga2   = 2*ga1;


F    = [];
G    = [];
H    = [];
feas = 0;

x23 = [x2;x3];
xs  = [x1;x23];
psi = x23.^a2;
psi = prod(psi);
x1s = x1.^2;
xi  = psi - x1s;
phi = psi./xi;
xm1 = 1./xs;

% first check if point is feasible. If not, stop immediately:
if any(xi < 0) || any(any(x23<0))
    feas = -1;
    return;
end

% Function value:
if compute(1) > 0
    F = - sum(log(xi)) - sum(sum(a1.*log(x23)));
end

% Gradient:
if compute(2) > 0 || compute(3) > 0
    tmp2    = [phi;phi;phi];
    tmp2    = (ga2.*tmp2 + 1 - ga1).*xm1;
    gt      = -tmp2;
    gt(1,:) = 2*x1 ./ xi;
    G       = gt(:);
end

% Hessian:
if compute(3) > 0

	% compute the 6 different elements in the 
	% 3x3 symmetric matrix one element at a time
	% Only need to compute the six different elements in here
	% and then output them as a vector as below. I.e. do not 
	% gather them into a block diagonal matrix.


    % element (2,3):
    Dphi  = (-x1s.*phi)./xi;
    Dphi3 = 2*((1-alpha').*Dphi).*xm1(3,:);
    el32  = -2*alpha'.*Dphi3.*xm1(2,:);

    % elements (2,2) and (3,3):
    Dphij = ([Dphi;Dphi;Dphi].*ga2).*xm1;    
    dd    = xm1.*(tmp2 - ga2.*Dphij);
    el22  = dd(2,:);
    el33  = dd(3,:);

    % elements (3,1) and (2,1)
    tmp  = -(4*x1.*phi./xi);
    col1 = [tmp;tmp].*ga1(2:3,:)./x23;	
    el21 = col1(1,:);
    el31 = col1(2,:);

    % element (1,1):
    el11 = (2./xi).*(1+2*x1s./xi);

    % output as vector in the following way:
    H = [el11;el21;el31;el21;el22;el32;el31;el32;el33];
    H = H(:);

end