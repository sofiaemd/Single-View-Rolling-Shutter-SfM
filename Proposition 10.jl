using Oscar
using LinearAlgebra

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
    C = vcat([0],[sum([c[j][i]*x^i for i in 1:d]) for j in 1:2])
    return C
end

#Taking the image of a line

C = makeC()
r = [1, 0, -x]

lineCurve = cross(r, cayleyTransform()*(q + cross(Delta, C)))

coeffs = vcat([coeff(lineCurve[1],i) for i in 1:2*delta+d+1],[coeff(lineCurve[2],i) for i in 0:2*delta+d+1])

# Evaluating and constructing the ideal

evalRandom = [rand(-20:20) for i in 1:3*delta+2*d+3]

evalCoeffs = [evaluate(p, evalRandom) for p in coeffs]

I = ideal(S, [coeffs[i] - evalCoeffs[i] for i in 1:4*delta+2*d+3])

print("The dimension of I is", dim(I))
print("The degree of I is", degree(I))

# Calculating its Jacobian in the same point as before

J = jacobian_matrix(coeffs)
Jdet = det(map(e -> evaluate(e, evalRandom), J[1:3*delta+2*d+3,1:3*delta+2*d+3]))


