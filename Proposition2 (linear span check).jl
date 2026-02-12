using AbstractAlgebra
using LinearAlgebra

R,x = polynomial_ring(QQ, :x)

function A(delta)
    if delta == 0
        return [0,0,0]
    else
        return [sum([rand(-20:20)*x^(i) for i=1:delta]) for j = 1:3]
    end
end

function cayleyTransform(b)
    # cayley transform without the normalization
    M=matrix(R,[[1+b[1]^2-b[2]^2-b[3]^2, 2*(b[1]*b[2]-b[3]), 2*(b[1]*b[3]+b[2])],
    [2*(b[1]*b[2]+b[3]), 1-b[1]^2+b[2]^2-b[3]^2, 2*(b[2]*b[3]-b[1])],
    [2*(b[1]*b[3]-b[2]),2*(b[2]*b[3]+b[1]),1-b[1]^2-b[2]^2+b[3]^2]])
    return M
end

function C(d)
    if d == 0
        return [0,0,0]
    else
        return [sum([rand(-20:20)*x^(i) for i=1:d]) for j = 1:3]
    end
end

function randomPlucker()
    points = rand(-20:20,6)
    return vcat(points[1:3]-points[4:6], cross(points[1:3],points[4:6]))
end

function lineCurve(d, delta, l)
    P = [cayleyTransform(A(delta)), C(d)]
    L = randomPlucker()
    r = [1,0,-x]
    u = cross(r,P[1]*(L[1:3]+cross(L[4:6],P[2])))
    coeffs1 = [coeff(u[1],i) for i = 1:(d+2*delta+1)]
    coeffs2 = [coeff(u[2],i) for i = 0:(d+2*delta+1)]
    coeffs = vcat(coeffs1,coeffs2)
    if l >= 2
    for i in 2:l
        L = randomPlucker()
        r = [1,0,-x]
        u = cross(r,P[1]*(L[1:3]+cross(L[4:6],P[2])))
        coeffs1 = [coeff(u[1],i) for i = 1:(d+2*delta+1)]
        coeffs2 = [coeff(u[2],i) for i = 0:(d+2*delta+1)]
        coeffs = vcat(coeffs,coeffs1,coeffs2)
    end
    end
    return coeffs
end

function isLinearSpanFull(d, delta, l)
    M = matrix([lineCurve(d, delta, l) for i in 1:l*(2*d+4*delta+3)])
    if det(M) != 0
        return true
    else 
        return false
    end
end

for d in 0:5
    for delta in 0:5
        for l in 0:6
            print(isLinearSpanFull(d,delta,l))
        end
    end
end