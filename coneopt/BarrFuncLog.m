function [F,G,H,feas] = BarrFuncLog(x1,x2,x3,K,want)

F     = [];
G     = [];
H     = [];
feas  = 0;

x2m1  = 1./x2;
x3m1  = 1./x3;
x2dx3 = x2./x3;
tmp1  = log(x2dx3);
psi   = x1 - x2.*tmp1;
psim1 = 1./psi;

% Feasibility:
if any(psi < 0) || any(any(x2<0)) || any(any(x3<0))
    feas = -1;
    return;
end

% Function value:
if want(1) > 0
    F = sum( -log(psi) - tmp1 );
end

if want(2) > 0
    tmp2  = [ -psim1;...
        psim1.*(1+tmp1)-x2m1;...
        -x3m1.*(x2.*psim1 + 1)];
    G = tmp2(:);
end

if want(3) > 0
    psi2  = psi.*psi;
    psim2 = 1./psi2;
    x2m2  = x2m1.*x2m1;
    x3m1  = 1./x3;
    x3m2  = x3m1.*x3m1;
    tmp3  = (1+tmp1);

    el11  = psim2;
    el21  = -tmp3.*psim2;
    el31  = x2dx3.*psim2;
    el22  = psim1.*x2m1 + x2m2 + psim2.*tmp3.^2;
    el32  = -x3m1.*psim1 - x2dx3.*psim2.*tmp3;
    el33  = x3m2.*(1+x2.*psim1 + psim2.*(x2.^2));

    H  = [el11;el21;el31;el21;el22;el32;el31;el32;el33];
    H  = H(:);

end