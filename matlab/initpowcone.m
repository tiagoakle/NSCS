function K = initpowcone(K,pars)

tmp         = [K.powc.alph';1-K.powc.alph'];
K.powc.a1   = flipud(tmp);
K.powc.a2   = 2*tmp;
K.powc.ga1  = [ones(1,K.npow);tmp];
K.powc.ga2  = 2*K.powc.ga1;
K.powc.a1d  = tmp;
K.powc.a23  = K.powc.a1d';
K.powc.a23  = K.powc.a23(:);
K.powc.apad = tmp.^tmp;

% "param-2 pow cone:"
K.powc.a1g2  = zeros(2,K.npow);
K.powc.ga1g2 = [ones(1,K.npow);zeros(2,K.npow)];

%glineur pow cone:
K.powc.a1g  = [1-2*K.powc.alph';2*K.powc.alph'-1];
K.powc.a1gm = max(K.powc.a1g, 0);
K.powc.ga1g = [ones(1,K.npow);K.powc.a1gm];

% "parabola" from 2 to 3, pow cone
K.powc.a1g3  = 4*[(1/2-K.powc.alph').^2;...
                   (K.powc.alph'-1/2).^2];
K.powc.ga1g3 = [ones(1,K.npow);K.powc.a1g3];

%nesterov pow cone:
K.powc.a1n  = ones(2,K.npow);
K.powc.ga1n = ones(3,K.npow);

K.powc.idx.x2x3 = [K.powc.idx.x2';K.powc.idx.x3'];

K.powc.fpownest = 0;
if pars.fpownest
    K.powc.fpownest = 1;
end
K.powc.fpowglin = 0;
if pars.fpowglin
    K.powc.fpowglin = 1;
end
K.powc.fpowglin2 = 0;
if pars.fpowglin2
    K.powc.fpowglin2 = 1;
end
K.powc.fpowglin3 = 0;
if pars.fpowglin3
    K.powc.fpowglin3 = 1;
end
K.powc.fpowhild = 0;
if pars.fpowhild
    K.powc.fpowhild = 1;
end
