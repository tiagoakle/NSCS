#Stores the problem definition 
type OpProblem
    m::Int
    n::Int
    A::Array{Float64,2}
    b::Array{Float64,1}
    c::Array{Float64,1}
    nu::Float64
end

#Constructor that initializes the problem definition
#and sets the complexity value
function OpProblem( m::Int,
                    n::Int,
                    A::Array{Float64,2},
                    b::Array{Float64,1},
                    c::Array{Float64,1})
    #TODO: Set nu to n+1, change this for LPs
    OpProblem(m,n,A,b,c,float64(n+1)) 
end

#This type stores the algorithm parameters
type OpPars
    maxIter::Int
    rho::Float64
    sigma::Float64
end

#Constructor that sets default parameters
function OpPars(prob::OpProblem)    
    #10 max iters, rho = nu+sqrt(nu)
    rho = prob.nu+sqrt(prob.nu)
    sigma = prob.nu/rho
    OpPars(10,rho,sigma)
end

#Stores the state of the solution
type OpState
    #Present point
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
    t::Float64
    k::Float64

    #Residuals
    d::Array{Float64,1} #dual residual
    p::Array{Float64,1} #Primal residual
    g::Float64          #Gap residual

    #Working vectors 
    rhs::Array{Float64,1} #To store the rhs as we build it 
    dir::Array{Float64,1} #To store the directions [dy;dx;dt]
    wk1::Array{Float64,1} #Used int the iterative refinement of the search directions
    ds::Array{Float64,1}  #To store the direction ds
    dk::Float64           #To store the direction for kappa

end

#Initialize the state with a given starting vector
function OpState(prob,
                 x::Array{Float64,1},
                 y::Array{Float64,1},
                 s::Array{Float64,1},
                 t::Float64,
                 k::Float64)

    #Allocate the residuals
    p = zeros(prob.m)
    d = zeros(prob.n)
    g = 0.0
    #Allocate the working space
    rhs = zeros(prob.m+prob.n+1)
    dir = zeros(prob.m+prob.n+1)
    wk1 = zeros(prob.m+prob.n+1)
    ds  = zeros(prob.n)
    dk = 0.0
    #TODO: verify belonging of x,s to the cones
    OpState(x,y,s,t,k,d,p,g,rhs,dir,wk1,ds,dk)
end

#Initialize the state object from the problem definition
#Sets the initial x,y,s,t,k to e,0,g(e),t,k
function OpState(prob::OpProblem)
    #TODO e should come from the cones 
    x = ones(prob.n)
    y = zeros(prob.m)
    s = ones(prob.n)
    t = 1.0
    k = 1.0
    OpState(prob,x,y,s,t,k)
end
