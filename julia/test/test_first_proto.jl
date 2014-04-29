#Tests the first prototype
include("../src/first_proto.jl")

using NSCS
#Generates a random matrix and tests the sovlver
m = 10
n = 100
A   = randn(m,n) #Generate the linear constraint matrix
b   = A*rand(n)     #Generate a rhs
c   = A'*randn(m) + rand(n) #Generate a bounded objective

#Form the problem structure
prob = NSCS.OpProblem(m,n,A,b,c)
NSCS.NSCSConeSolver(prob)
