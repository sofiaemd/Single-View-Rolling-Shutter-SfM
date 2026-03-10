using Oscar
using LinearAlgebra
using HomotopyContinuation


# Functions for camera

function cayleyTransform(Svars,x,d,delta)
    a = [[Svars[i+3*j+2*d-4] for j in 1:delta] for i in 1:3]
    A = [sum([a[j][i]*x^i for i in 1:delta]) for j in 1:3]
    return matrix([
        1+A[1]^2-A[2]^2-A[3]^2  2*(A[1]*A[2]-A[3])    2*(A[1]*A[3]+A[2]);
        2*(A[1]*A[2]+A[3])      1-A[1]^2+A[2]^2-A[3]^2 2*(A[2]*A[3]-A[1]);
        2*(A[1]*A[3]-A[2])      2*(A[2]*A[3]+A[1])    1-A[1]^2-A[2]^2+A[3]^2
    ])
end


function makeC(Svars,x,d,chart)
    if chart == 1
        c = [vcat(1, Svars[2:d]),vcat(Svars[1], Svars[d+1:2*d-1])]
    elseif chart == 2
        c = [vcat(Svars[1], Svars[2:d]),vcat(1, Svars[d+1:2*d-1])]
    end
    C = vcat([0],[sum([c[j][i]*x^(d-i+1) for i in 1:d]) for j in 1:2])
    return C
end

# Function converting pols to use homotopy continuation

function hcPols(l,y)
    return [sum((Float64(numerator(c)) / Float64(denominator(c))) * prod(y[i]^e[i] for i in 1:length(e)) for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))) for f in l]
end

function toFloat(l)
    return [Float64(numerator(c)) / Float64(denominator(c)) for c in l]
end

function toComplex(list)
    return [ComplexF64(numerator(c)) / ComplexF64(denominator(c)) for c in list]
end

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

function lineCurve(Svars,x,d,delta,q,Delta,chart)
    C = makeC(Svars,x,d,chart)
    r = [1, 0, -x]
    curve = cross(r, cayleyTransform(Svars,x,d,delta)*(q + cross(Delta, C)))
    return vcat([-1*coeff(curve[1],i) for i in 1:2*delta+d+1],[-1*coeff(curve[2],i) for i in 1:2*delta+d+1])
end

function nrSols(d, delta, coeffs, evalRandom, evalCoeffs)
    @var y[1:2*d+3*delta+3]
    @var p[1:4*delta+2*d+2]

    coeffsHc = hcPols(coeffs,y)
    coeffsWParam = [coeffsHc[i]-p[i] for i in 1:4*delta+2*d+2]

    F = System(coeffsWParam, variables = y, parameters = p)

    subF = System(coeffsWParam[1:2*d+3*delta+3], variables = y, parameters = p[1:2*d+3*delta+4])
    sols = monodromy_solve(subF, toFloat(evalRandom), toFloat(evalCoeffs[1:2*d+3*delta+4]))
    solsF = countSols(F, solutions(sols),toFloat(evalCoeffs))
    return solsF
end

function checkProp10(d,delta)
    print("***************  ","Case d = ",d," and delta = ", delta,"  ***************", "\n")
    S, Svars = polynomial_ring(QQ, 3*delta+2*d+3)
    R, x = polynomial_ring(S, :x)
    Delta = [1,Svars[2*d+3*delta],Svars[2*d+3*delta+1]]
    q = [-Svars[2*d+3*delta]*Svars[2*d+3*delta+2]-Svars[2*d+3*delta+1]*Svars[2*d+3*delta+3],Svars[2*d+3*delta+2],Svars[2*d+3*delta+3]]

    #Obtaining the polynoials for the curve coefficients and evaluating for a specific camera and line
    coeffs = lineCurve(Svars,x,d,delta,q,Delta,1)

    evalRandom = vcat([1],[rand(vcat(-20:-1, 1:20)) for i in 1:3*delta+2*d+2])
    evalCoeffs = [AbstractAlgebra.evaluate(p, evalRandom) for p in coeffs]
    print("The parameters of the chosen camera and line are: ", evalRandom,"\n")
    print("The image of the chosen point through the picture-taking map: ", evalCoeffs,"\n")

    #Finding the solutions in the first chart
    print("Checking if there are any other solutions in the chart b_d = 1","\n")
    s = nrSols(d, delta, coeffs, evalRandom, evalCoeffs)
    if length(s) == 1
        if s[1]==toComplex(evalRandom)
            print("The chosen point is the only solution to the system")
        else
            print("The solution found does not match the original point")
        end
    else
        print(length(s), "solutions have been found")
    end
    print("\n ----------------------- \n")

    #Finding the soltions in the second chart
    print("Checking if there are any other solutions in the chart c_d = 1","\n")
    coeffs = lineCurve(Svars,x,d,delta,q,Delta,2)
    s = nrSols(d, delta, coeffs, evalRandom, evalCoeffs)
    if length(s) == 1
        if s[1]==toComplex(evalRandom)
            print("The chosen point is the only solution to the system")
        else
            print("The solution found does not match the original point")
        end
    else
        print(length(s), "solutions have been found")
    end
    print("\n ----------------------- \n")
        
    #jacobian check
    J = jacobian_matrix(coeffs)
    Jval = map(e -> AbstractAlgebra.evaluate(e,evalRandom),J)
    if rank(Jval) == 3*delta+2*d+3
        print("The map is full rank at the chosen point")
    else
        print("The map is not full rank at the chosen point")
    end

    print("\n\n\n")
end


# First case takes less than 10 minutes, second case takes around 4 hours and third case about a day
checkProp10(1,2)
checkProp10(2,2)
checkProp10(3,2)

