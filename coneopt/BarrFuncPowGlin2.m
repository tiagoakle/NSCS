function [F,G,H,feas] = BarrFuncPowGlin2(x1,x2,x3,K,want)

% -log(x^(2a)*y^(2-2a) - z^2)

F    = [];
G    = [];
H    = [];
feas = 0;

x23 = [x2;x3];
xs  = [x1;x23];

psi = x23.^K.powc.a2;
psi = prod(psi);
x1s = x1.^2;
xi  = psi - x1s;
phi = psi./xi;
xm1 = 1./xs;

% feasibility:
if any(xi < 0) || any(any(x23<0))
    feas = -1;
    return;
end

% Function value:
if want(1) > 0
    F = - sum(log(xi)) - sum(sum(K.powc.a1g2.*log(x23)));
end

% Gradient:
if want(2) > 0 || want(3) > 0
    tmp2    = [phi;phi;phi];
    tmp2    = (K.powc.ga2.*tmp2 + K.powc.ga1g2).*xm1;
    gt      = -tmp2;
    gt(1,:) = 2*x1 ./ xi;
    G       = gt(:);
end


% Hessian:
if want(3) > 0

    % element (2,3):
    Dphi  = (-x1s.*phi)./xi;
    Dphi3 = 2*((1-K.powc.alph').*Dphi).*xm1(3,:);
    el32  = -2*K.powc.alph'.*Dphi3.*xm1(2,:);

    % elements (1,1), (2,2) and (3,3) (DIAGONAL):
    Dphij = ([Dphi;Dphi;Dphi].*K.powc.ga2).*xm1;
    dd    = xm1.*(tmp2 - K.powc.ga2.*Dphij);
    el22  = dd(2,:);
    el33  = dd(3,:);

    % first row and column:
    tmp  = -(4*x1.*phi./xi);
    col1 = [tmp;tmp].*K.powc.ga1(2:3,:)./x23;
    el21 = col1(1,:);
    el31 = col1(2,:);

    % element (1,1):
    el11 = (2./xi).*(1+2*x1s./xi);

    % gather it all:
    H = [el11;el21;el31;el21;el22;el32;el31;el32;el33];
    H = H(:);

end


