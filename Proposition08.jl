using Oscar #to define the polynomial system
using LinearAlgebra #to do cross product
using HomotopyContinuation #to solve the system using monodromy


# Functions for camera

function cayleyTransform(Svars,x,delta,d)
    a = [[Svars[3*d-1 + j + (delta)*(i-1)] for j in 1:delta] for i in 1:3]
    A = [sum([a[j][i]*x^i for i in 1:delta]) for j in 1:3]
    return matrix([
        1+A[1]^2-A[2]^2-A[3]^2  2*(A[1]*A[2]-A[3])    2*(A[1]*A[3]+A[2]);
        2*(A[1]*A[2]+A[3])      1-A[1]^2+A[2]^2-A[3]^2 2*(A[2]*A[3]-A[1]);
        2*(A[1]*A[3]-A[2])      2*(A[2]*A[3]+A[1])    1-A[1]^2-A[2]^2+A[3]^2
    ])
end


function makeC(Svars,x,chart,d)
    if chart == 1
        c = [[i + j == 2 ? 1 : Svars[i+(j-1)*3-1] for j in 1:d] for i in 1:3]
    elseif chart == 2
        c = [[i + j == 2 ? Svars[1] : i==2 && j == 1 ? 1 : Svars[i+(j-1)*3-1] for j in 1:d] for i in 1:3]
    elseif chart == 3
        c = [[i + j == 2 ? Svars[2] : i==3 && j == 1 ? 1 : Svars[i+(j-1)*3-1] for j in 1:d] for i in 1:3]
    end
    C = [sum([c[j][i]*x^(d-i+1) for i in 1:d]) for j in 1:3]
    return C
end

# Function converting polinomials and number lists to use homotopy continuation

function hcPols(list,y)
    return [sum((Float64(numerator(c)) / Float64(denominator(c))) * prod(y[i]^e[i] for i in 1:length(e)) for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))) for f in list]
end

function toFloat(list)
    return [Float64(numerator(c)) / Float64(denominator(c)) for c in list]
end

function toComplex(list)
    return [ComplexF64(numerator(c)) / ComplexF64(denominator(c)) for c in list]
end

#Function that returns the polynomials that define the coefficients of the image curve given a camera and a line

function lineCurve(Svars,x,line,chart,d,delta)
    Delta = [line[1],line[2],-(line[1])*(line[3])-(line[2])*(line[4])]
    q = [line[3],line[4],1] 
    C = makeC(Svars,x,chart,d)
    r = [1, 0, -x]

    curve = cross(r, cayleyTransform(Svars,x,delta,d)*(q + cross(Delta, C)))

    coeffs = vcat([-1*coeff(curve[1],i) for i in 1:2*delta+d+1],[-1*coeff(curve[2],i) for i in 1:2*delta+d+1])
    return coeffs
end

# Function that goes through the solutions of the square polynomial subsystem and verifies 
# which if them satisfy the whole polynomial system. 
# It returns a list of the solutions that satisfy the whole system

function countSols(F,sols, parameters)
    l = []
    tol = 1e-8
    for s in sols
    a = maximum(abs.(HomotopyContinuation.evaluate(F, s,parameters)))
    if a < tol
        push!(l,s)
    end
    end
    return l
end




function checkprop8(d,delta)


print("***************  ","Case d = ",d," and delta = ", delta,"  ***************", "\n")

S, Svars = polynomial_ring(QQ, 3*delta+3*d+7)
R, x = polynomial_ring(S, :x)

# list of variables corresponding to the two lines
L = [Svars[3*d+3*delta+i*4:3*d+3*delta+4*i+3] for i in 0:1]

# Choosing a generic camera and lines. 
# We choose a_d, b_d and c_d to be 1 and all the other coefficients to be a random non zero integer between -30 and 30
randomEval = vcat([1,1],rand(vcat(-30:-1, 1:30),3*delta+3*d+5))
print("The parameters of the chosen camera and two lines are: ", randomEval,"\n")

# We construct the polynomials defining the coefficients for the picture-taking map and evaluate them on the chosen point
coeffs = [lineCurve(Svars,x,L[i],1,d,delta) for i in 1:2]
evalCoeffs = [[AbstractAlgebra.evaluate(p, randomEval) for p in coeffs[j]] for j in 1:2]
print("The image of the chosen point through the picture-taking map: ", evalCoeffs,"\n")


# Variables to convert the polynomial system so it can be solved using HomotopyContinuation
@var y[1:3*delta+3*d+7]
@var p[1:(4*delta+2*d+2)*2]

# We will check that the two curves define a unique camera by searching for solutions in the charts where
# a_d=1, b_d=1 and c_d=1, this will be enough to ensure uniqueness in the whole space. Our previous choice of 
# camera coefficients ensures that the curve coefficients are the same in the three charts. 
params = toFloat(vcat(evalCoeffs[1],evalCoeffs[2]))

for chart in 1:3
# We construct the polynomial system that defines the picture-taking map
coeffs = [lineCurve(Svars,x,L[i],chart,d,delta) for i in 1:2]
coeffsHc = [hcPols(e,y) for e in coeffs]
coeffWParam = vcat([coeffsHc[1][i]-p[i] for i in 1:4*delta+2*d+2],[coeffsHc[2][i]-p[4*delta+2*d+2+i] for i in 1:4*delta+2*d+2])
F = System(coeffWParam, variables = y, parameters = p)

if chart == 1
    print("Checking if there are any other solutions in the chart a_d = 1","\n")
elseif chart == 2
    print("Checking if there are any other solutions in the chart b_d = 1","\n")
else
    print("Checking if there are any other solutions in the chart c_d = 1","\n")
end

# Using monodromy to find all the solutions to a square subsystem of F
subCoeffWParam = vcat([coeffsHc[1][i]-p[i] for i in 1:2*delta+d+3],[coeffsHc[2][i]-p[2*delta+d+3+i] for i in 1:delta+2*d+4])
subF = System(subCoeffWParam, variables = y, parameters = p[1:3*d+3*delta+7])
paramsSub = toFloat(vcat(evalCoeffs[1][1:2*delta+d+3],evalCoeffs[2][1:delta+2*d+4]))
sols = monodromy_solve(subF, toFloat(randomEval), paramsSub)

# Checking that the only solution that satisfies the whole system is the chosen point
s = countSols(F, solutions(sols), params)
if length(s) == 1
    if s[1]==toComplex(randomEval)
        print("The chosen point is the only solution to the system")
    else
        print("The solution found does not match the original point")
    end
else
    print("Multiple solutions have been found")
end

print("\n ----------------------- \n")

end

# Calculating the Jacobian in the chosen point
J = jacobian_matrix(reduce(vcat,coeffs))
Jval = map(e -> AbstractAlgebra.evaluate(e,randomEval),J)
if rank(Jval) == 3*delta+3*d+7
    print("The map is full rank at the chosen point")
else
    print("The map is not full rank at the chosen point")
end

print("\n\n\n")
end

checkprop8(1,1)
checkprop8(1,2)
checkprop8(2,1)
