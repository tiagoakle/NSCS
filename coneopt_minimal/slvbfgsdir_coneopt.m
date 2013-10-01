function [d,v] = slvbfgsdir_coneopt(CF,v,rs,upd,pars)

% get the Choleschy factors:
L   = CF.L;
N   = CF.N;

% get the current right hand side:
r1  = rs.r1;
r2  = rs.r2;
r3  = rs.r3;
r4  = rs.r4;
r5  = rs.r5;

% update the BFGS data holders:
idx1              = 2*v.j-1;
idx2              = idx1+1;
v.Q1(:,idx1:idx2) = [upd.v,upd.w];
v.sc1(idx1:idx2)  = [upd.sgnv;upd.sgnw];
v.Q2(:,idx1:idx2) = pars.A*v.Q1(:,idx1:idx2);
if pars.permuteM
    tmp                = pars.Sp'*v.Q2(:,idx1:idx2);
    tmp                = N\(N'\tmp);
    v.MQ2(:,idx1:idx2) = pars.Sp*tmp;
else
    v.MQ2(:,idx1:idx2) = N\(N'\(v.Q2(:,idx1:idx2)));
end
Q3               = v.Q2(:,1:idx2);
MQ3              = v.MQ2(:,1:idx2);
sc3              = v.sc1(1:idx2);
MS               = (diag(sc3) + Q3'*MQ3);

% solve sequence:=======
sc1r  = repmat(v.sc1,1,2);
r6    = r2 + r5;
tmp   = [-r6,pars.c];
% tmp = B*tmp
tmp   = L\(L'\tmp) + v.Q1*(sc1r.*(v.Q1'*tmp));
    
tmp   = [r1,pars.b] + pars.A*tmp;
% this is tmp = M*tmp, where M ~= (A*H^{-1}*A')^{-1}
tmp   = tmp - Q3*( MS\(MQ3'*tmp) );
if pars.permuteM
    tmp = pars.Sp'*tmp;
end
tmp   = N\(N'\tmp);
if pars.permuteM
    tmp = pars.Sp*tmp;
end
% --------------------
r9    = tmp(:,1);
r10   = tmp(:,2);
tmp   = [r6,-pars.c] + pars.A'*tmp;
% tmp = B*tmp
tmp   = L\(L'\tmp) + v.Q1*(sc1r.*(v.Q1'*tmp));
r11   = tmp(:,1);
r12   = tmp(:,2);
% for ease of notation (and copying later):
c     = pars.c;
b     = pars.b;
A     = pars.A;
tau   = v.tau;
kappa = v.kappa;
% ------------------
t1    = r4 + tau*(c'*r11 - b'*r9 + r3);
t2    = tau*(b'*r10 - c'*r12) + kappa;
d{2}  = t1/t2;
d{5}  = (r4 - kappa*d{2})/tau;
d{3}  = r9  + r10*d{2};
d{1}  = r11 + r12*d{2};
d{4}  = c*d{2} - r2 - A'*d{3};




