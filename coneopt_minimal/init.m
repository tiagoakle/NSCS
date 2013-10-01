function [v,K,pars,R] = init(varargin)

% A, b, c, v0, K, pars

if nargin < 5
    error('coneopt requires at least 5 input: A,b,c,v0,K.')
end

if nargin > 5
    pars = varargin{6};
else
    pars = [];
end

% get input:
pars.A = varargin{1};
pars.b = varargin{2};
pars.c = varargin{3};
v0     = varargin{4};
K      = varargin{5};

% set the defaults:
pars   = setdefaults(pars);

% initialize the cone:
K      = initcone(K,pars);

% check initial points:
if ~isfield(v0,'x');
    error('v0 must contain at least field x)');
end

v.x = v0.x;

v.F = BarrFuncP(v.x,K,[1,1,1]);
if v.F{4} < 0
    error('initial point x0 is NOT in primal cone!');
end

if isfield(v0,'tau');
    v.tau = v0.tau;
else
    v.tau = 1;
end

if isfield(v0,'kappa');
    v.kappa = v0.kappa;
else
    v.kappa = 1/v.tau;
end

if isfield(v0,'s');
    v.s = v0.s;
else
    %v.s = -v.mu*v.F{2};

    % choosing s so that
    % s + (x'*s+tau*kappa)/(nu+1) * grad F(x) = 0
    % this is just a set of linear equations
    % solve them for s -- everything else is known
    % the system is easy to solve because
    % the system matrix is (I+v*x'). Use
    % Woodbury formula:
    % ---------------------------------------
    vt  = v.F{2}/(K.nu+1);
    qt  = -v.tau*v.kappa*v.F{2}/(K.nu+1);
    v.s = qt-(v.x'*qt)*vt/(1+v.x'*vt);
end

% first scaling point:
v.u = v.x;

if isfield(v0,'mu');
    v.mu = v0.mu;
else
    v.mu = (v.x'*v.s + v.tau*v.kappa)/(K.nu+1);
end
v.mu0 = v.mu;

% this is not allowed when requiring that we not use 
% the dual barrier:
%FD = BarrFuncD(v.s,K,[1,-1,-1]);
%if FD{4} < 0
%    error('initial point s0 is NOT in dual cone!');
%end

if isfield(v0,'y');
    v.y = v0.y;
else
    v.y = ones(pars.m,1);
end


% stopping constants:
pars.relstopP  = max(1,norm([pars.A,pars.b],'inf'));
pars.relstopD  = max(1,norm([pars.A',speye(pars.n),-pars.c],'inf'));
pars.relstopG  = max(1,norm([-pars.c',pars.b',1],'inf'));

% duality gap:
% v.dgap = v.x'*v.s + v.tau*v.kappa;

% allocate space for the "reordered" Hessian into 
% blockdiag matrix of 1x1 (pos) and 3x3 (other 3-dim cones):
% ...

% precision level:
if isfield(pars,'preclvl')
    
    if (round(pars.preclvl)==pars.preclvl) && pars.preclvl > 0
    
        rhoPs = [ 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-6 ];
        rhoDs = [ 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-6 ];
        rhoAs = [ 1e-1, 1e-3, 1e-4, 1e-5, 1e-8,  1e-6 ];
        rhoGs = [ 1e-1, 1e-3, 1e-5, 1e-6, 1e-8,  1e-6 ];
        rhoIs = [ 1e-2, 1e-3, 1e-5, 1e-6, 1e-8,  1e-8 ];
        rhoMs = [ 1e-4, 1e-6, 1e-8, 1e-9, 1e-10, 1e-8 ];
        oit   = [ 300,  300,  300,  500,  1000,  500  ];
        iit   = [ 300,  300,  300,  500,  1000   500  ];
        
        pars.rhoP = rhoPs(pars.preclvl);
        pars.rhoD = rhoDs(pars.preclvl);
        pars.rhoA = rhoAs(pars.preclvl);
        pars.rhoG = rhoGs(pars.preclvl);
        pars.rhoI = rhoIs(pars.preclvl);
        pars.rhoM = rhoMs(pars.preclvl);
    
        pars.outermaxit = oit(pars.preclvl);
        pars.innermaxit = iit(pars.preclvl);
        
    end
end


if pars.permuteM
    if ~issparse(pars.A)
        fprintf('init: warning: A NOT sparse! \n');
        fprintf('      setting pars.permuteM = 0 \n');
        pars.permuteM = 0;
    else
        % generate the permutation of the matrix
        % M corresponding to the sparsity pattern of M:
        L           = chol(v.F{3});   
        KT          = (L')\(pars.A');  
        M           = KT'*KT;     
        [jnk,jnk,S] = chol(M);
        pars.Sp     = S;
    end
end

% other init:
R.v0              = v;
R.pars            = pars;
R.K               = K;
R.dat.nchols      = 0;
R.dat.nls         = 0;
R.dat.ncentsteps  = 0; % total
v.i               = 0;
v.j               = 0;
v.jt              = 0;
R.dat.nkktforms   = 0;
R.dat.nkktsolves  = 0;
R.status          = 'not done';
R.stop            = false;
v.lam             = 99;
v.a               = 0;

if abs(pars.echo) > 1
    if pars.trace < 1
        error('if |echo|>1, must have trace > 0');
    end
end

% BFGS data holders:
v.Q1  = zeros(pars.n,2*pars.cnbfgsstps);
v.Q2  = zeros(pars.m,2*pars.cnbfgsstps);
v.MQ2 = zeros(pars.m,2*pars.cnbfgsstps);
v.sc1 = zeros(2*pars.cnbfgsstps,1);

% do everything the first time:
v.k    = 1;
[v,R]  = getresiduals(v,pars,R);
R      = checkstopcrit(v,pars,R);
v.dgap = v.mu * (K.nu+1);
v.feas = 0;
R      = gettrace(v,pars,R);
R.trace.centstps(v.k) = 0;
R.trace.nbfgsstepst(v.k) = 0;
v.dtauaff = 0;
v.dkappaaff = 0;
v.tauprev = 1;
v.kappaprev = 1;
v.reachedinnermaxit = 0;

% printing:
printfunc('initial',v,K,pars,R);


tic;
R.dat.t1 = toc;


warning('off','MATLAB:nearlySingularMatrix');


