using Oscar
using LinearAlgebra
using HomotopyContinuation

#Defining the parameters and variables 

d = 1
delta = 2

S, Svars = polynomial_ring(QQ, 3*delta+2*d+3)
R, x = polynomial_ring(S, :x)

c = [[i + j == 2 ? 1 : Svars[d*(i-1)+j-1] for j in 1:d] for i in 1:2]
a = [[Svars[i+3*j+2*d-4] for j in 1:delta] for i in 1:3]
Delta = [1,Svars[2*d+3*delta],Svars[2*d+3*delta+1]]
q = [-Svars[2*d+3*delta]*Svars[2*d+3*delta+2]-Svars[2*d+3*delta+1]*Svars[2*d+3*delta+3],Svars[2*d+3*delta+2],Svars[2*d+3*delta+3]]

# Functions for camera

function cayleyTransform()
    A = [sum([a[j][i]*x^i for i in 1:delta]) for j in 1:3]
    return matrix([
        1+A[1]^2-A[2]^2-A[3]^2  2*(A[1]*A[2]-A[3])    2*(A[1]*A[3]+A[2]);
        2*(A[1]*A[2]+A[3])      1-A[1]^2+A[2]^2-A[3]^2 2*(A[2]*A[3]-A[1]);
        2*(A[1]*A[3]-A[2])      2*(A[2]*A[3]+A[1])    1-A[1]^2-A[2]^2+A[3]^2
    ])
end


function makeC()
    C = vcat([0],[sum([c[j][i]*x^(d-i+1) for i in 1:d]) for j in 1:2])
    return C
end

function makeC1()
    c0 = [[i + j == 2+d ? 1 : Svars[d*(i-1)+j] for j in 1:d] for i in 1:2]
    C = vcat([0],[sum([c0[j][i]*x^(i) for i in 1:d]) for j in 1:2])
    return C
end

# Function converting pols to use homotopy continuation

function hcPols(l)
    return [sum((Float64(numerator(c)) / Float64(denominator(c))) * prod(y[i]^e[i] for i in 1:length(e)) for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))) for f in l]
end

function toFloat(l)
    return [Float64(numerator(c)) / Float64(denominator(c)) for c in l]
end

#Taking the image of a line

C = makeC()
r = [1, 0, -x]

lineCurve = cross(r, cayleyTransform()*(q + cross(Delta, C)))

coeffs = vcat([coeff(lineCurve[1],i) for i in 1:2*delta+d+1],[coeff(lineCurve[2],i) for i in 0:2*delta+d+1])

# Evaluating and constructing the ideal

evalRandom = [rand(vcat(-20:-1, 1:20)) for i in 1:3*delta+2*d+3]

evalCoeffs = [AbstractAlgebra.evaluate(p, evalRandom) for p in coeffs]

# Doing Homotopy continuation to find solutions

@var y[1:2*d+3*delta+3]
@var p[1:4*delta+2*d+3]

coeffsHc = hcPols(coeffs)
coeffsWParam = [coeffsHc[i]*p[1]-p[i]*coeffsHc[1] for i in 2:4*delta+2*d+3]


F = System(coeffsWParam, variables = y, parameters = p)
sol = HomotopyContinuation.solve(F, target_parameters=toFloat(evalCoeffs))
print(sol)

# Checking on other chart

lineCurve1 = cross(r, cayleyTransform()*(q + cross(Delta, makeC1())))
coeffs1 = vcat([coeff(lineCurve1[1],i) for i in 1:2*delta+d+1],[coeff(lineCurve1[2],i) for i in 0:2*delta+d+1])
coeffsHc1 = hcPols(coeffs1)
coeffsWParam1 = [coeffsHc1[i]*p[1]-p[i]*coeffsHc1[1] for i in 2:4*delta+2*d+3]


F1 = System(coeffsWParam1, variables = y, parameters = p)
sol1 = HomotopyContinuation.solve(F1, target_parameters=toFloat(evalCoeffs))
print(sol1)

# TO DO: For minimal problems, check that the solutions is both charts are the same


# Calculating its Jacobian in the same point as before

J = jacobian_matrix(coeffs)
Jdet = det(map(e -> AbstractAlgebra.evaluate(e, evalRandom), J[1:3*delta+2*d+3,1:3*delta+2*d+3]))


