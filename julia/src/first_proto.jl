
module NSCS
include("types.jl")

#Calculates the largest step size that keeps x,t,k feasible and mu positive
function maxStep(state::OpState,problem::OpProblem)
    #Shorthand 
    m = problem.m
    n = problem.n
    
    x = state.x
    y = state.y
    s = state.s
    t = state.t
    k = state.k

    ds = state.ds
    dy = sub(state.dir,1:n) #Access these by reference
    dx = sub(state.dir,n+1:n+m)
    dt = state.dir[m+n+1]
   
    #Step size calculation
    #Find the maximum step that will keep x positive and mu positive
    dbx   = minimum(-x[dx.<0]./dx[dx.<0]) #Take the indexes that reduce the value of x
    dbt   = -t/dt
    if(dbt[1]>0) dbx = minimum((dbt[1],dbx)) end 
    alpha = dbx

    #dbk   = -k/dk
    #if(dbk[1]>0) dbx = minimum((dbk[1],dbx)) end 
    
    #Calculate the step size that would keep mu positive
    pm_a  = dx'*ds+dt*dk
    pm_a  = pm_a[1]
    pm_b  = dx'*s + dt*k+ds'*x+dk*t
    pm_b  = pm_b[1]
    if(pm_b.^2-4*pm_a>0)
        root1 = (-pm_b+sqrt(pm_b.^2-4*pm_a))/(2*pm_a)
        root2 = (-pm_b-sqrt(pm_b.^2-4*pm_a))/(2*pm_a)
        if root1 < 0 root1 = Inf end
        if root2 < 0 root2 = Inf end 
     
        #Rules to decide the maximum step length
        alpha = minimum((dbx,root1,root2))
    end
    if alpha == Inf alpha = 1 end
    alpha = 0.99*alpha   
end


#Generates a random matrix and tests the sovlver
m = 10
n = 100
A   = randn(m,n) #Generate the linear constraint matrix
b   = A*rand(n)     #Generate a rhs
c   = A'*randn(m) + rand(n) #Generate a bounded objective

prob = OpProblem(m,n,A,b,c)
pars = OpPars(prob)

#Inirial value of the variables
x   = rand(n)
s   = rand(n)
y   = randn(m)
t   = rand()
k   = rand()

mu  = x'*s+t*k
mu  = mu/nu 
mu  = mu[1]


#Main iteration
for iter = 1:maxIter

    H   = diagm(1./(x.*x))
    h   = t^-2
    
    #Present value of the residuals
    p   = A*x-t.*b
    d   = -A'*y-s+t*c
    g   = b'*y-c'*x-k
    
    #Assemble the matrix 
    G = [zeros(m,m)   A     -b ;
              -A'     mu*H   c ;
               b'    -c'  mu*h ] 
    
    #Assemble the rhs
    rhs= [-p;
          -d-s+sigma*mu./x;
          -g-k+sigma*mu/t]
    
    #Solve!
    (L,U,P) = lu(G) #Factorize
    di    = L\rhs[P]
    dz    = U\di

    #Iterative refinement 
    for iref = 1:10
        di    = rhs-G*dz
        di    = di[P]
        di    = L\di
        dz    = U\di + dz
    end

    #Extract the step directions
    dy    = dz[1:m]
    dx    = dz[m+1:m+n]
    dt    = dz[m+n+1]
    ds    = -s+sigma*mu./x - mu*H*dx
    dk    = -k+sigma*mu/t - mu*h*dt
    
    #REMOVE:Evaluate the error
    #pe    = A*dx-dt*b+p
    #de    = -A'*dy+dt*c-ds+d
    #ge    = b'*dy-c'*dx-dk+g
    #dxe   = mu*H*dx+ds+s-sigma*mu./x
    #dte   = mu*h*dt+dk+k-sigma*mu/t
    #println("np $(norm(pe)) nd $(norm(de)) ng $ge nxe $(norm(dxe)) nte $(abs(dte))")

    alpha = maxStep(problem,state)
    #Take the step and validate that it is in the cone and mu is positive

    y = y+alpha*dy
    x = x+alpha*dx
    t = t+alpha*dt
    s = s+alpha*ds
    k = k+alpha*dk
    
    #Update mu 
    mu = (x'*s+t.*k)/nu
    mu = mu[1]

    #Validate the step did not cause an infeasibility
    if minimum(x) < 0 println("Infeasible x") end
    if t[1] < 0 println("Infeasible t") end
    #if k[1] < 0 println("Infeasible k") end
    if mu < 1.e-15
        break
    end
    if mu < 0 println("Mu negative") 
        break 
    end 
    println("iter $iter, alpha $alpha, dbx $dbx, dbt $dbt, mu $mu")
end
end
