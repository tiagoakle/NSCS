#Generates a random matrix and tests the sovlver
m = 10
n = 100
A   = randn(m,n) #Generate the linear constraint matrix
b   = A*rand(n)     #Generate a rhs
c   = A'*randn(m) + rand(n) #Generate a bounded objective

#Algorithm parameters 
nu      = n+1
rho     = nu+sqrt(nu)
sigma   = nu/rho
maxIter = 1

#Inirial value of the variables
x   = rand(n)
s   = rand(n)
y   = randn(m)
t   = rand(1)
k   = rand(1)

#Main iteration
    mu  = x'*s+t.*k
    mu  = mu/(n+1) 

    H   = diagm(x.*x)
    h   = t.^2
    
    #Present value of the residuals
    p   = A*x-t.*b
    d   = -A'*y-s+t.*c
    g   = b'*y-c'*x-k
    
    #Assemble the matrix 
    G = [zeros(m,m)   A  -b ;
              -A'     H   c ;
               b'    -c'  h ] 
    
    #Assemble the rhs
    rhs= [-p;
          -d-s+sigma*mu./x;
          -g-k+sigma*mu/t]
    
    #Solve!
    (L,U) = lu(G) #Factorize
    di    = L\rhs 
    dz    = U\di
    #Iterative refinement 
    di    = L\(+rhs-G*dz)
    dz    = U\di + dz
    #Extract the step directions
    dy    = dz[1:m]
    dx    = dz[m+1:m+n]
    dt    = dz[m+n+1]
    ds    = -s+sigma*mu./x - H*dx
    dk    = -k+sigma*mu./t - h*dt
    
    #Step size calculation
    #Find the maximum step that will keep x positive and mu positive
    ix    = 0
    dbx   = minimum(x[dx.<0]) #Take the indexes that reduce the value of x
    pm_a  = dx'*ds+dt*dk
    pm_b  = dx'*s + dt.*k+ds'*x+dk.*t
    root1 = (-pm_b+sqrt(pm_b.^2-4*pm_a))/(2*pm_a)
    root2 = (-pm_b-sqrt(pm_b.^2-4*pm_a))/(2*pm_a)
    println(root2)


