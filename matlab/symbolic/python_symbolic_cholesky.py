#USes sympy to compute the symbolic cholesky factors of H
import sympy
from sympy import *
x = Symbol('x')
y = Symbol('y',positive=True)
z = Symbol('z',positive=True) 

f = -sympy.log(z*sympy.log(y/z)-x)-sympy.log(y)-sympy.log(z);

h11 = f.diff(x,x)
h21 = f.diff(y,x)
h22 = f.diff(y,y)
h31 = f.diff(z,x)
h32 = f.diff(z,y)
h33 = f.diff(z,z)

H = sympy.Matrix([[h11,h21,h31],[h21,h22,h32],[h31,h32,h33]]);
# h11 h21 h31 
# h21 h22 h23
# h31 h32 h33

c11 = sympy.sqrt(h11);
c21 = h21/c11;
c31 = h31/c11;

#Calculate the remainder 2by2 matrix
r22 = h22-c21*c21
r32 = h32-c21*c31
r33 = h33-c31*c31

#Calculate the next column of the cholesky
c22 = sympy.sqrt(r22)
c32 = r32/c22

#calculate the remainder 
r33 = r33-c32*c32
c33 = sympy.sqrt(r33)

#Test
C = sympy.Matrix([[c11,0,0],[c21,c22,0],[c31,c32,c33]])
D = C*C.transpose()-H
D.simplify()


r = Symbol('r')
l = Symbol('l')
C.simplify()
#C[0,0] = 1/r
C = C.subs(-x+z*sympy.log(y/z),r)
C = C.subs(sympy.log(y/z),l)

E = C.subs(l,sympy.log(y/z))
E = E.subs(r,-x+z*sympy.log(y/z))
D = E*E.transpose()-H
D.simplify()

SiC = Matrix([[1/r,0,0],[-z/(r*y),sympy.sqrt((r+z)/(r*y**2)),0],[(-l+1)/r, -1/sympy.sqrt(r**2+r*z), sympy.sqrt((r+2*z)/(z**2*(r+z)))]])


S = sympy.Matrix([
[                1/r,  0,   0],
[    -z/(r*y), (y**(-2) + z/(r*y**2))**(1/2),                                                                                                           0],
[-(l - 1)/(r), -1/(r*y*(y**(-2) + z/(r*y**2))**(1/2)), (z**(-2) + 1/(r*z) - (-l + 1)*(l - 1)/r**2 - (l - 1)**2/r**2 - 1/(r**2*y**2*(y**(-2) + z/(r*y**2))))**(1/2)]])

E = S.subs(l,log(y/z))
E = E.subs(r,-x+z*log(y/z))
D = E*E.transpose()-H
D.simplify()

T  = sympy.Matrix( [[1/r,0,0], [-z/(r*y),y*(1+z/r)**(1/2),0],[-(l-1)/r,-1/(r*(1+z/r)**(1/2)), 1/r*(r/z*(1+r/z)-1/(y**2*(1+z/r)))]] )

E = T.subs(l,log(y/z))
E = E.subs(r,-x+z*log(y/z))
D = E*E.transpose()-H
D.simplify()

T  = sympy.Matrix( [[1/r,0,0], [-z/(r*y),1/y*(1+z/r)**(1/2),0],[-(l-1)/r,-1/(r*(1+z/r)**(1/2)), (1/z**2+1/(r*z)-1/(r**2*(1+z/r)))**(1/2)]] )

E = T.subs(l,log(y/z))
E = E.subs(r,-x+z*log(y/z))
D = E*E.transpose()-H
D.simplify()

T  = sympy.Matrix( [[1/r,0,0], [c21,c22,0],[c31,c32,c33]])

E = T.subs(l,log(y/z))
E = E.subs(r,-x+z*log(y/z))
D = E*E.transpose()-H
D.simplify()

