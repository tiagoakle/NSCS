function K = initcone(K,pars)

% barrier parameter:
expcoef = 3;
logcoef = 3;
poscoef = 1;
if pars.fpownest
    powcoef = 4*K.npow;
elseif pars.fpowglin
    aa = K.powc.alph;
    powcoef = sum(max(3-2*aa , 1+2*aa));
elseif pars.fpowglin2
    powcoef = 2*K.npow;
elseif pars.fpowglin3
    aa = K.powc.alph;
    powcoef = sum(2+4*(1/2-aa).^2);
elseif pars.fpowhild
    ps  = 1./K.powc.alph;
    qs  = 1./(1 - 1./ps);
    mpq = max(ps,qs);
    powcoef = sum( 3*mpq./(mpq+1) );
else
    powcoef = 3*K.npow;
end
K.nu = poscoef*K.npos + powcoef + ...
    expcoef*K.nexp + logcoef*K.nlog;

% generate indices:
ndiffcones = 0;
if K.npos > 0
    K.posc.idx.x = 1:K.npos;
    ndiffcones   = ndiffcones + 1;
else
    K.posc.idx.x = [];
end

if K.npow > 0
    tmp           = 1:K.npow;
    K.powc.idx.x1 = K.npos + tmp;
    K.powc.idx.x2 = K.powc.idx.x1(end) + tmp;
    K.powc.idx.x3 = K.powc.idx.x2(end) + tmp;
    ndiffcones    = ndiffcones + 1;
else
    K.powc.idx.x1 = [];
    K.powc.idx.x2 = [];
    K.powc.idx.x3 = [];
end

if K.nexp > 0
    tmp           = 1:K.nexp;
    K.expc.idx.x1 = K.npos + 3*K.npow  + tmp;
    K.expc.idx.x2 = K.expc.idx.x1(end) + tmp;
    K.expc.idx.x3 = K.expc.idx.x2(end) + tmp;
    ndiffcones    = ndiffcones + 1;
else
    K.expc.idx.x1 = [];
    K.expc.idx.x2 = [];
    K.expc.idx.x3 = [];
end

if K.nlog > 0
    tmp           = 1:K.nlog;
    K.logc.idx.x1 = K.npos + 3*K.npow + 3*K.nexp + tmp;
    K.logc.idx.x2 = K.logc.idx.x1(end) + tmp;
    K.logc.idx.x3 = K.logc.idx.x2(end) + tmp;
    ndiffcones    = ndiffcones + 1;
else
    K.logc.idx.x1 = [];
    K.logc.idx.x2 = [];
    K.logc.idx.x3 = [];
end

% subscripts for accumulation of function val:
K.allsubsF = ones(ndiffcones,2);

% Subscripts for accumulation of gradient and hessian of barrfunc

% pos cone
K.posc.idx.subsG = (1:K.npos)';
K.posc.idx.subsH = [K.posc.idx.subsG,K.posc.idx.subsG];

% pow cone
tmp = [K.powc.idx.x1;K.powc.idx.x2;K.powc.idx.x3];
K.powc.idx.subsG = tmp(:);
if K.npow > 0
    sb1        = [tmp;tmp;tmp];
    sb2        = [tmp(1,:);tmp(1,:);tmp(1,:);...
        tmp(2,:);tmp(2,:);tmp(2,:);...
        tmp(3,:);tmp(3,:);tmp(3,:)];
else
    sb1 = []; sb2 = [];
end
K.powc.idx.subsH = [sb1(:),sb2(:)];

% exp cone:
tmp = [K.expc.idx.x1;K.expc.idx.x2;K.expc.idx.x3];
K.expc.idx.subsG = tmp(:);
if K.nexp > 0
    sb1        = [tmp;tmp;tmp];
    sb2        = [tmp(1,:);tmp(1,:);tmp(1,:);...
        tmp(2,:);tmp(2,:);tmp(2,:);...
        tmp(3,:);tmp(3,:);tmp(3,:)];
else
    sb1 = []; sb2 = [];
end
K.expc.idx.subsH = [sb1(:),sb2(:)];

% log cone:
tmp = [K.logc.idx.x1;K.logc.idx.x2;K.logc.idx.x3];
K.logc.idx.subsG = tmp(:);
if K.nlog > 0
    sb1        = [tmp;tmp;tmp];
    sb2        = [tmp(1,:);tmp(1,:);tmp(1,:);...
        tmp(2,:);tmp(2,:);tmp(2,:);...
        tmp(3,:);tmp(3,:);tmp(3,:)];
else
    sb1 = []; sb2 = [];
end
K.logc.idx.subsH = [sb1(:),sb2(:)];


% gather all:
K.allsubsG = [K.posc.idx.subsG;K.powc.idx.subsG;...
    K.expc.idx.subsG;K.logc.idx.subsG];
K.allsubsH = [K.posc.idx.subsH;K.powc.idx.subsH;...
    K.expc.idx.subsH;K.logc.idx.subsH];

% generate LINEAR indices for the HESSIAN:
n = pars.n; n = [n,n];
K.Hidxl = sub2ind(n,K.allsubsH(:,1),K.allsubsH(:,2));
K.posc.idx.subsHl = sub2ind(n,K.posc.idx.subsH(:,1),K.posc.idx.subsH(:,2));
K.powc.idx.subsHl = sub2ind(n,K.powc.idx.subsH(:,1),K.powc.idx.subsH(:,2));
K.expc.idx.subsHl = sub2ind(n,K.expc.idx.subsH(:,1),K.expc.idx.subsH(:,2));
K.logc.idx.subsHl = sub2ind(n,K.logc.idx.subsH(:,1),K.logc.idx.subsH(:,2));

% in blocks of 9 in each column (hessians):
K.powc.idx.subsHlblk = reshape(K.powc.idx.subsHl,9,K.npow);
K.expc.idx.subsHlblk = reshape(K.expc.idx.subsHl,9,K.nexp);
K.logc.idx.subsHlblk = reshape(K.logc.idx.subsHl,9,K.nlog);

% in blocks of 3 in each column (gradients and variables):
K.powc.idx.subsGlblk = reshape(K.powc.idx.subsG,3,K.npow);
K.expc.idx.subsGlblk = reshape(K.expc.idx.subsG,3,K.nexp);
K.logc.idx.subsGlblk = reshape(K.logc.idx.subsG,3,K.nlog);

% and gathered:
K.powexpc.idx.subsGlblk = [K.powc.idx.subsGlblk,K.expc.idx.subsGlblk,...
    K.logc.idx.subsGlblk];
K.powexpc.idx.subsHlblk = [K.powc.idx.subsHlblk,K.expc.idx.subsHlblk,...
    K.logc.idx.subsHlblk];

% extra for power cone (pre-compute some stuff with alphas):
if K.npow > 0
    K = initpowcone(K,pars);
end
% extra for exp cone:
if K.nexp > 0
    K = initexpcone(K,pars);
end




