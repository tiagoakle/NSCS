%Symbolic calculation of the tensor product
% y = \nabla^3f(x)[H^{-1}a,H^{-1}a]

%Define the symbols
syms x y z a1 a2 a3 real
f = -log(z*log(y/z)-x)-log(y)-log(z);
a = [a1;a2;a3];
H = hessian(f);
Hi = inv(H);
%Calculate the three parts of the tensor with the jacobian
J = jacobian(H(:));
T1 = reshape(J(:,1),3,3);
T2 = reshape(J(:,2),3,3);
T3 = reshape(J(:,3),3,3);

y1 = a'*Hi*T1*Hi*a;
y2 = a'*Hi*T2*Hi*a;
y3 = a'*Hi*T3*Hi*a;

y1 = simplify(y1,2000)
y2 = simplify(y2,2000)
y3 = simplify(y3,2000)
