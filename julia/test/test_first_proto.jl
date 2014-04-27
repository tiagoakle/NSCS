#Generates a problem and uses first proto to solve it

#Tests for the first proto file
include("../src/types.jl")

#Build a problem file
m = 10
n = 100
A   = randn(m,n) #Generate the linear constraint matrix
b   = A*rand(n)     #Generate a rhs
c   = A'*randn(m) + rand(n) #Generate a bounded objective

#Test that the complexity is being calculated correctly
prob = OpProblem(m,n,A,b,c)
@assert prob.nu == n+1

#Test that the default parameters are set correctly
pars = OpPars(prob)
@assert pars.rho == prob.nu+sqrt(prob.nu)
@assert pars.sigma == prob.nu/pars.rho

#Test the OpState Constructors
state = OpState(prob)
#state = OpState(prob,randn(n),rand(m),rand(n),1.0,1.0)
@assert state.p == zeros(m)
@assert state.d == zeros(n)
@assert state.g == 0.0
@assert state.rhs == zeros(m+n+1)
@assert state.dir == zeros(m+n+1)
@assert state.wk1 == zeros(m+n+1)
@assert state.ds == zeros(n)
@assert state.dk == 0.0

#Test the OpState Constructors
state = OpState(prob,randn(n),rand(m),rand(n),1.0,1.0)
@assert state.p == zeros(m)
@assert state.d == zeros(n)
@assert state.g == 0.0
@assert state.rhs == zeros(m+n+1)
@assert state.dir == zeros(m+n+1)
@assert state.wk1 == zeros(m+n+1)
@assert state.ds == zeros(n)
@assert state.dk == 0.0

