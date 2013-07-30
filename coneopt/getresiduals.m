function [v,R] = getresiduals(v,pars,R)

bty  = pars.b'*v.y;
ctx  = pars.c'*v.x;

v.rA = abs( ctx - bty )/( v.tau + abs(bty) );

v.rP = (pars.A*v.x - pars.b*v.tau);
v.rD = (-pars.A'*v.y - v.s + pars.c*v.tau);
v.rG = (-ctx + bty - v.kappa);

v.rPrel = norm( v.rP, 'inf')/pars.relstopP;
v.rDrel = norm( v.rD, 'inf')/pars.relstopD;
v.rGrel = norm( v.rG, 'inf')/pars.relstopG;
v.rArel = v.rA;

%v.rH = (v.kappa*v.tau - v.mu);
%v.rC = (v.s + v.mu*v.F{2});