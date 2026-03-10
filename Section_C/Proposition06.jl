using Oscar

l = 6
R, x = polynomial_ring(QQ, 5+4*l)

# This function returns the polynomials for the coefficients we will use
function coefficients(a,b,c,D,q)
    Z1 = -q[1]+D[2]*a[1]-D[1]*b[1]
    Zlast = D[3]*b[2]-D[2]*c[2]
    Y0 = q[2]
    Y1 = D[3]*a[1]-D[1]*c[1]
    Ylast = D[3]*a[2]-D[1]*c[2]
    return [Z1,Zlast,Y0,Y1,Ylast]
end

#Variables in the chart c_d=1
a = x[1:2]
b = x[3:4]
c = [x[5],1]
lines = [x[6+4*i:9+4*i] for i in 0:l-1]

#Coefficients in the chart c_d=1
eqs = [coefficients(a,b,c,[L[1],L[2],-L[1]*L[3]+L[2]*L[4]],[L[3],L[4],1]) for L in lines]

# Point X where we will evaluate
randomEval = [rand(vcat(-20:-1, 1:20)) for i in 1:5+4*l]

#Ideal corresponding to the points that yield the same coefficients as X
I = ideal(reduce(vcat,[[p - evaluate(p,randomEval) for p in e] for e in eqs]))
dim(I)
degree(I)

#Checking again in the chart c_d=0
c = [x[5],0]
eqs = [coefficients(a,b,c,[L[1],L[2],-L[1]*L[3]+L[2]*L[4]],[L[3],L[4],1]) for L in lines]
I = ideal(reduce(vcat,[[p - evaluate(p,randomEval) for p in e] for e in eqs]))
dim(I)
degree(I)