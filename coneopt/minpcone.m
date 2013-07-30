function [x,f,R] = minpcone(A,b,p,varargin)

% min_x  || x ||_p
%  s.t.    Ax = b
%
% see page 16 in paper (Skajaa, Ye)


if nargin > 3
    pars = varargin{1};
end

[M,N]   = size(A);
b       = b(:);
M2      = length(b);

if M ~= M2
    error('A and b do not match');
end
if p < 1
    error('p must be >= 1');
end
if M > N
    error('M must be <= N');
end

% scale A and b:
maxAb = max(max([A,b]));
A     = A/maxAb;
b     = b/maxAb;

pars.n = 3*N;
pars.m = M + 1 + N - 1;

% build cone:
K.npos = 0;
K.npow = N;
K.nexp = 0;
K.nlog = 0;

% power-cone specifics:
K.powc.alph = (1/p)*ones(K.npow,1);

% build A:
AA                    = zeros(pars.m,pars.n);
AA(1:M,1:N)           = A;
AA(M+1,N+1:2*N)       = 1;
AA(M+1,2*N+1)         = -1;
AA(M+2:end,2*N+1:end) = ...
    full(spdiags([ones(N-1,1),-ones(N-1,1)],[0 1],N-1,N));

% it does NOT appear to be faster to make A sparse:
AA = sparse(AA);
%

% build b:
bb      = zeros(pars.m,1);
bb(1:M) = b;

% build c:
cc        = zeros(pars.n,1);
cc(2*N+1) = 1;

% starting point:
e     = ones(N,1);
% x0    = A\b+0.1*ones(N,1);
% tmp   = abs(x0);
% tmp   = ((N*e).^K.powc.alph).*tmp;
% uk    = max(tmp);
% y0    = uk./e;
% u0    = uk.*(N*e);

x0 = zeros(N,1);
y0 = e;
u0 = e;

v0.x  = [x0;y0;u0];

% call to coneopt:
R = coneopt(AA,bb,cc,v0,K,pars);

if strcmp(R.status,'optimal')
    x = R.sol.xopt(1:N);
    f = norm( x,p );
else
    x = [];
    f = inf;
end

