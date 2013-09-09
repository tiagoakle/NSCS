function [F,G,H,feas] = BarrFuncPos(x,K,want)

F    = [];
G    = [];
H    = [];
feas = 0;

% feasibility positive cone:
if any(x<0)
    feas = -1;
    return;
end

if want(1) > 0
    F = -sum(log(x));
end
if want(2) > 0
    G = -1./x;
end
if want(3) > 0
    H = 1./x.^2;
end