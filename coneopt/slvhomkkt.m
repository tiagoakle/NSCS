function [d,CF] = slvhomkkt(H,mu,A,b,c,tau,kappa,r1,r2,r3,r4,r5,pars)

% using structure, "efficiently" solve the system:
%
% A*dx - b*dtau           = r1
% -A'*dy - ds + c*dtau    = r2
% -c'*dx + b'*dy - dkappa = r3
% tau*dkappa + kappa*dtau = r4
% (mu*H)*dx + ds          = r5
%
% x, tau, y, s, kappa

[L,d{6}] = chol(H);    % Choleschy L'*L = H
if d{6} > 0            % if Choleschy failed,
    CF = [];           % rounding errors are preventing
    return;            % further progress, so exit
end                    %

L     = sqrt(mu)*L;    % L is then Choleschy of mu*H
if nargout > 1
    CF.L = L;
end

KT    = (L')\(A');     % these two lines
M     = KT'*KT;        % time critical
r6    = r2 + r5;       %

tmp   = L'\([-r6,c]);
tmp   = [r1,b] + KT'*tmp;

% factorization
if pars.permuteM
    M = pars.Sp'*(M*pars.Sp);
end
[N,d{6}] = chol(M);
if d{6} > 0
    if pars.cholinc
        N = cholinc(M,'inf');
    else
        return;
    end
end
if pars.permuteM
    tmp = pars.Sp'*tmp;
    tmp = N\(N'\tmp);
    tmp = pars.Sp*tmp;
else
    tmp   = N\(N'\tmp);
end

% OLD VERSION
% to preserve as much sparsity in chol(M) as possible:
%if pars.permuteM
%
%    Mp       = pars.Sp'*(M*pars.Sp);
%    [N,d{6}] = chol(Mp);   % also time critical
%    if d{6} > 0            % if Choleschy failed,
%        return;            % rounding errors are preventing
%    end                    % further progress, so exit
% 
%     tmp = pars.Sp'*tmp;
%     tmp = N\(N'\tmp);
%     tmp = pars.Sp*tmp;
% 
% else
% 
%     [N,d{6}] = chol(M);
%     if d{6} > 0            % if Choleschy failed,
%         return;            % rounding errors are preventing
%     end                    % further progress, so exit
%     tmp   = N\(N'\tmp);
% 
% end

if nargout > 1
    CF.N = N;
end

r9    = tmp(:,1);
r10   = tmp(:,2);
tmp   = [r6,-c] + A'*tmp;
tmp   = L\(L'\tmp);
r11   = tmp(:,1);
r12   = tmp(:,2);

t1    = r4 + tau*(c'*r11 - b'*r9 + r3);
t2    = tau*(b'*r10 - c'*r12) + kappa;

d{2}  = t1/t2;
d{5}  = (r4 - kappa*d{2})/tau;
d{3}  = r9  + r10*d{2};
d{1}  = r11 + r12*d{2};
d{4}  = c*d{2} - r2 - A'*d{3};




