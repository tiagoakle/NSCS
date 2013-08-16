function upd = bfgsupd(s,y,By)

% check the curvature condition:
ys  = s'*y;
if ys <= 0
    fprintf('%9.5e \n',ys);
    error('error in quasi-Newton update');
end
rho = 1/ys;

% BFGS update:
Q     = orth([s,By]);
Qts   = Q'*s;
QtBy  = Q'*By;
Q2    = Qts*QtBy';
Z     = rho*(rho * y'*By + 1)*(Qts*Qts') - rho*(Q2 + Q2');

[V,L] = eig(Z);
evA   = Q*V;

cv    = L(1,1);
cw    = L(2,2);
v     = evA(:,1);
w     = evA(:,2);

upd.sgnv = sign(cv);
upd.v    = sqrt(abs(cv))*v;

upd.sgnw = sign(cw);
upd.w    = sqrt(abs(cw))*w;


% if cv < 0
%     upd.sv   = '-';
%     upd.sgnv = -1;
% else
%     upd.sv   = '+';
%     upd.sgnv = 1;
% end
% cv     = abs(cv);
% upd.v  = sqrt(cv)*v;
% 
% if cw < 0
%     upd.sw   = '-';
%     upd.sgnw = -1;
% else
%     upd.sw   = '+';
%     upd.sgnw = 1;
% end
% cw     = abs(cw);
% upd.w  = sqrt(cw)*w;

