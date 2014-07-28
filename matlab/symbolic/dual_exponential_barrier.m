%Calculate the hessian, and gradient for the barrier of the dual cone
clear all
x = sym('x','positive');
y = sym('y','positive');
z = sym('z','real');

r = sym('r','positive');
l = sym('l','real');

%The real barrier for the dual cone is 
%f = -log(-x*log(-y/x)-x+z)-log(y)-log(-x), however since the symbolic
%system only allows for postitie constraints we use 
%the substitution x for -x
f = -log(x*log(-y/-x)+x+z)-log(y)-log(x)
H = hessian(f);
g = gradient(f);

%Incorporate the assumptions of positivity
Hs = subs(H,l,log(y/x))
Hs = subs(Hs,r,x*l+x+z)
%Hs = subs(Hs,-x*l-x+z,-r);

h11 = Hs(1,1);
h21 = Hs(2,1);
h31 = Hs(3,1);
h22 = Hs(2,2);
h32 = Hs(3,2);
h33 = Hs(3,3);

c11 = sqrt(h11);
c21 = h21/c11;
c31 = h31/c11;

%Calculate the remainder 2by2 matrix
r22 = h22-c21*c21;
r32 = h32-c21*c31;
r33 = h33-c31*c31;

%Calculate the next column of the cholesky
c22 = sqrt(r22);
c32 = r32/c22;

%calculate the remainder 
r33 = r33-c32*c32;
c33 = sqrt(r33);

%Test 
Cs = [[c11 0    0 ];...
     [c21 c22  0 ];...
     [c31 c32 c33]];
Cs = simplify(Cs)

E= Hs-Cs*Cs.';
E = subs(E,log(y/x),l)
E = subs(E,x*l+x+z,r)
simplify(E)

%The substitution is not working 