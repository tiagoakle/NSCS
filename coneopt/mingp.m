function [x,f,R] = mingp(A,d,B,k,varargin)

if nargin > 4
    pars = varargin{1};
end

% some input checks here:
% 
% ---------

% P is the number of posynomial constraints
% M is the number of monomial constraints
% N is the number of variables
% NJ(j) is the number of monomials in posynomial j

P     = size(A,1);  
N     = size(A{1},2);
M     = size(B,1);    
NJ    = zeros(P,1); 
for j = 1:P         
    NJ(j) = size(A{j},1);
end

if isempty(B) || isempty(k)
    B = zeros(0,N);
end

Nh     = sum(NJ);
pars.n = 3*Nh + 2*N + 1;
pars.m = 2*Nh + M + P;

% build cone:
K.npos = 2*N+1;
K.npow = 0;
K.nexp = Nh;
K.nlog = 0;
K      = getbarrpar(K);

% build A, b, c:
%c
cc    = sparse([],[],[],pars.n,1,1);
cc(1) = 1;

%b
bb                    = sparse([],[],[],pars.m,1,pars.m);
bb(2:(P+Nh))          = 1;
h                     = log(k);
bb((P+Nh+1):(P+Nh+M)) = -h;
c                     = log(cell2mat(d));
bb((P+Nh+M+1):pars.m) = -c;

%A
nnzAA = 2*Nh+2*M*N+2*Nh*N+Nh+1;
AA = sparse([],[],[],pars.m,pars.n,nnzAA);
AA(1,1) = -1;
tmp = cell(P,1);
for j = 1:P
    tmp{j} = ones(1,NJ(j));
end
AA(1:P,1+2*N+Nh+1:1+2*N+2*Nh)           = blkdiag(tmp{:});
AA(P+1:P+Nh,1+2*N+2*Nh+1:1+2*N+3*Nh)    = speye(Nh);
AA(P+Nh+1:P+Nh+M,2:N+1)                 = B;
AA(P+Nh+1:P+Nh+M,N+1+1:2*N+1)           = -B;
AA(P+Nh+M+1:P+Nh+M+Nh,2*N+1+1:2*N+1+Nh) = -speye(Nh);
AA(P+Nh+M+1:P+Nh+M+NJ(1),2:2*N+1)       = [A{1},-A{1}];

% THIS IS SLOW!!
for j = 2:P
    idxnow = P+Nh+M+sum(NJ(1:j-1))+1:P+Nh+M+sum(NJ(1:j));
    AA(idxnow,2:2*N+1) = [A{j},-A{j}];
end


% starting point:
t00 = 1;
up0 = ones(N,1);
um0 = ones(N,1);
w00 = -ones(Nh,1);
v00 = ones(Nh,1);  
y00 = 0.5*ones(Nh,1);

v0.x  = [t00;up0;um0;w00;v00;y00];

% call to coneopt:
R = coneopt(AA,bb,cc,v0,K,pars);

% extract solution here:
if strcmp(R.status,'optimal')
    up  = R.sol.xopt(2:N+1);
    um  = R.sol.xopt(N+2:2*N+1);
    u   = up-um;
    x   = exp(u);
    f   = R.sol.xopt(1); % because cc=[1,0,0,...,0].
    % another way to get f:
    %xr  = repmat(x',NJ(1),1);
    %tmp = xr.^A{1}; 
    %tmp = prod(tmp,2);
    %f  = d{1}'*tmp;
else
    x = [];
    f = inf;
end

