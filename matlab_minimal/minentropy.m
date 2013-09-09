function [x,f,R] = minentropy(A,b,d,varargin)

% min_x   sum_j  d_j x_j log(x_j) 
% s.t.     Ax = b
%           x >= 0
%
% see page 18 in paper (Skajaa, Ye)

if nargin > 3
    pars = varargin{1};
end

[M,N]   = size(A);
b       = b(:);
M2      = length(b);

if M ~= M2
    error('A and b do not match');
end

pars.n = 3*N;
pars.m = N+M;

% build cone:
K.npos = 0;
K.npow = 0;
K.nexp = N;
K.nlog = 0;
K      = getbarrpar(K);

% build A:
nnzA    = nnz(A);
AA      = sparse([],[],[],pars.m,pars.n,nnzA+N);
tmp     = N*(N+M);
tmp     = tmp + 1:(N+M+1):(2*N*(N+M)-M);
AA(tmp) = 1;
tmp     = N+1:pars.m;
tmp2    = 2*N+1:pars.n;
AA(tmp,tmp2) = A;

% build b:
bb             = ones(pars.m,1);
bb(N+1:pars.m) = b;
bb             = sparse(bb);

% build c:
cc      = zeros(pars.n,1);
cc(1:N) = -d;
cc      = sparse(cc);

% starting point:
u0  = -ones(N,1);
v00 = ones(N,1);  
x0  = 0.5*ones(N,1);

v0.x  = [u0;v00;x0];

% special pars:
% it appears to be better not to permute
% ATLEAST FOR THE SPECIFIC PROBLEM FROM ERLING!
if ~isfield(pars,'permuteM')
    pars.permuteM = 0; 
end

% call to coneopt:
R = coneopt(AA,bb,cc,v0,K,pars);

% extract solution here:
if strcmp(R.status,'optimal')
    x = R.sol.xopt(2*N+1:pars.n);
    u = R.sol.xopt(1:N);
    f = -d'*u;
else
    x = [];
    f = inf;
end



