function [F,G,H,feas] = BarrFuncExp(x1,x2,x3,K,want)

F     = [];
G     = [];
H     = [];
feas  = 0;

logx2 = log(x2);
logx3 = log(x3);
x2m1  = 1./x2;
x3m1  = 1./x3;
tmp1  = logx2-logx3;
psi   = x3.*tmp1 - x1;
psim1 = 1./psi;
xi    = tmp1 - 1;

% Feasibility:
if any(psi < 0) || any(any(x2<0)) || any(any(x3<0))
    feas = -1;
    return;
end

% Function value:
if want(1) > 0
    F = sum( -log(psi) - logx2 - logx3 );
end

if want(2) > 0
    tmp2  = [ psim1;...
        -x2m1.*(x3.*psim1 + 1);...
        -xi.*psim1 - x3m1];
    G = tmp2(:);
end

if want(3) > 0
    psi2  = psi.*psi;
    psim2 = psim1.*psim1;
    x2m2  = x2m1.*x2m1;
    x3m1  = 1./x3;
    x3m2  = x3m1.*x3m1;

    el11  = psim2;
    el21  = -x3.*x2m1.*psim2;
    el31  = -xi.*psim2;
    el22  = psim2.*x2m2.*(x3.*psi + x3.*x3 + psi2);
    el32  = psim2.*x2m1.*(x3.*xi - psi);
    el33  = psim2.*(x3m1.*psi + xi.*xi + psi2.*x3m2);

    H  = [el11;el21;el31;el21;el22;el32;el31;el32;el33];
    H  = H(:);

end