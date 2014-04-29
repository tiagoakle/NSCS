module NSCS
include("types.jl")
using Debug
using Base.BLAS

#Calculates the largest step length that keeps x,t,k feasible and mu positive
function maxStep(state::OpState,problem::OpProblem)
    #Shorthand variables
    m = problem.m
    n = problem.n
    
    #Refferences to the present point
    x  = state.x
    t  = state.t
    s  = state.s
    k  = state.k
    
    #Refferences to the search direction
    ds = state.ds
    dy = sub(state.dir,1:m) 
    dx = sub(state.dir,m+1:n+m)
    dt = state.dir[m+n+1]
    dk = state.dk
   
    #Step size calculation

    #Find the maximum step that will keep x and s positive 
    dbx = Inf
    for j=1:n
        tmp = -x[j]/dx[j]
        if(dx[j]<0 && tmp < dbx)
            dbx = tmp
        end
        
        tmp = -s[j]/ds[j]
        if(ds[j]<0 && tmp < dbx)
            dbx = tmp
        end
    end 
    alpha = dbx
    
    #Find the step size that will also keep tau and kappa positive
    dbt   = -t[1]/dt[1]
   
    if(dbt>0) alpha = minimum((dbt,alpha)) end 

    dbk   = -k[1]/dk[1]
    if(dbk>0) alpha = minimum((dbk,alpha)) end 
    
    #Calculate the step size that would keep mu positive
    #Observe that x's will become x's+alpha(dx's+ds'x) + alpha^2 dx'ds 
    pm_a  = dot(n,pointer(state.dir)+sizeof(Float64)*m,1,ds,1)+ dt*dk
    pm_b  = dot(n,pointer(state.dir)+sizeof(Float64)*m,1,s,1) + dot(n,ds,1,x,1) + dt*k + dk*t
    
    println(pm_a)
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

#Print the header message
function PrintHeader(prob::OpProblem,pars::OpPars)
if(!pars.silent)
    println("=============================================")
    println("      Non-symmetric cone solver              ")
    println(" Pos: $(prob.n), linear constraints $(prob.m)")
    println(" Santiago Akle: SOL Stanford Univerity 2014  ")
    println(" 2014                                        ")
    println("=============================================")
end
end

#Calculates the initial point for the cone problem 
#sets e for all the symmetric cones, the "centered" 
#point for the exponential, sets y to zeros
function SetInitialPoint(problem::OpProblem)
    #TODO: add the initialization for the non symmetric cone
    #TODO: add the initialization for SOCP
    #Allocate the space for x,y
    x = ones(problem.n)
    y = zeros(problem.m)
    (x,y)
end

#Computes the initial dual s
function ComputeInitialS(prob::OpProblem,x::Array{Float64,1})
    #TODO: Initial s for nonsymmetric cones
    s = ones(prob.n)
    s
end

#Execute the solution of the problem with the default initial point 
# and default parameters
function NSCSConeSolver(prob::OpProblem)
    pars = OpPars(prob)
    #TODO: Validate the problem entries
    x,y = SetInitialPoint(prob)
    #Call the solver
    NSCSConeSolver(prob,x,y,pars)
end

#Exectute the solution of the problem with a prescribed initial primal x 
#and dual y. The dual s will be determined from the gradient of the barrier at x
#and default parameters
function NSCSConeSolver(prob::OpProblem,x::Array{Float64,1},y::Array{Float64,1},pars::OpPars)
     #Initialize the state with x,y,s,t=1.0,k=1.0
    s     = ComputeInitialS(prob,x)
    state = OpState(prob,y,x,1.0,s,1.0)
   
    #Shorthand for the problem definition
    m   = prob.m
    n   = prob.n
    nu  = prob.nu
    A   = prob.A
    b   = prob.b
    c   = prob.c
    
    rho = pars.rho
    sigma = pars.sigma
    maxIter = pars.maxIter
    
    #Calculate the initial mu
    state.mu = dot(n,state.x,1,state.s,1) + state.t*state.k
    state.mu = state.mu/nu

    #Print header
    PrintHeader(prob,pars)

    #Main iteration
    for iter = 1:maxIter
        
        #Evaluate the hessians 
        H   = diagm(1./(state.x.*state.x))
        h    = state.t.^-2
        #H    = diagm(state.s./state.x)
        #h    = state.k/state.t
        
        #Present value of the residuals
        #TODO: Investivgate if A,b,c here are copies or references
        state.p   = A*state.x-state.t.*b
        state.d   = -A'*state.y-state.s+state.t.*c
        state.g   = (b'*state.y-c'*state.x-state.k)[1]
        
        #Assemble the matrix 
        G = [zeros(m,m)   A                 -b ;
                  -A'     state.mu*H         c ;
                  b'    -c'         state.mu*h ] 
                 #-A'               H         c ;
                 # b'    -c'        h           ] 
        
        #Assemble the rhs
        state.rhs= [-state.p;
                    -state.d-state.s+sigma*state.mu./state.x;
                    -state.g-state.k+sigma*state.mu/state.t]
        
        #Solve!
        (L,U,P) = lu(G) #Factorize
        state.wk1     = L\state.rhs[P]
        state.dir     = U\state.wk1
        
        #Iterative refinement 
        for iref = 1:10
           state.wk1    = state.rhs-G*state.dir
           state.wk1    = state.wk1[P]
           state.wk1    = L\state.wk1
           state.dir    = U\state.wk1 + state.dir
        end

        #Calculate ds and dk
        state.ds    = -state.s+sigma*state.mu./state.x - state.mu*H*state.dir[m+1:m+n];
        state.dk    = (-state.k+sigma*state.mu./state.t - state.mu*h.*state.dir[m+n+1])[1]
        
        #Take the step and validate that it is in the cone and mu is positive
        alpha = maxStep(state,prob)
            
        #state.y = state.y+alpha*dy this would create a new vector for the add and multiply... :(
        axpy!(alpha,state.dir,1:m,state.y,1:m)     
        axpy!(alpha,state.dir,m+1:m+n,state.x,1:n) #axpy with format alpha, x, range, y , range
        state.t = state.t+alpha*state.dir[m+n+1]
        axpy!(alpha,state.ds,1:n,state.s,1:n)
        state.k = state.k+alpha*state.dk
        
        #Update mu 
        state.mu = dot(n,state.x,1,state.s,1) + state.t*state.k
        state.mu = state.mu/nu
        

        println("iter $iter, alpha $alpha, mu $(state.mu), min(x) $(minimum(x))")
        #Validate the step did not cause an infeasibility
        if minimum(state.x) < 0 println("Infeasible x, min(x) $(minimum(state.x))") end
        if state.mu < 0 println("Mu negative") 
            break 
        end 
        if state.t[1] < 0 println("Infeasible t") end
        #if k[1] < 0 println("Infeasible k") end
        if state.mu < 1.e-15 println("Mu reached 1.e-15")
            break
        end
        
    end #End main iteration
end #End function

end #Module NSCS
