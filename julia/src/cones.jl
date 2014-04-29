#Defines the functions for the different cones

#These functions are indicators, they take the value true if the point is feasible
function PosPrimalIn(x::Array{Float64})
end

function PosDualIn(x::Array{Float64})
end

function SOCPPrimalIn(x::Array{Float64})
end

function SOCPDualIn(x::Array{Float64})
end

function ExpPrimalIn(x::Array{Float64})
end

function ExpDualIn(x::Array{Float64})
end

#Set ident functions set the identity element for the 
#self-dual cones and the "centered" element for the non-self-dual cones
#for the cone of size n at the position x

#Set the n entries of x equal to one 
function SetIdentPos(n,x::Array{Float64,1})
    for j=1:n
        x[j] = 1.0
    end
end

#Set the n entries of the socp cone to e
function SetSOCPPos(n,x::Array{Float64,1})
    x[1] = 1.0
    for j=2:n
        x[j] = 0.0
    end
end

