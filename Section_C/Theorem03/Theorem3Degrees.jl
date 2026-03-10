using HomotopyContinuation

#The code presented here provides the numerical degree computation for the p=1 case. 
#Do not run the code with p>1. You can adapt d but be aware that it takes approx 8-9h for d=2.
d=1


p=1
@var x[1:2*(6*d+2)]
@var t0
V=x[1:6*d+2]
U=x[6*d+3:2*(6*d+2)]


#define the necessary functions,i.e. Cayley and the paramtrising polyomials
function cayleyNotNorm(a,b,c)
    return [1+a^2-b^2-c^2 2*(a*b-c) 2*(a*c+b); 2*(a*b+c) 1-a^2+b^2-c^2 2*(b*c-a); 2*(a*c-b) 2*(b*c+a) 1-a^2-b^2+c^2]
end

function alpha(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[i]*t^i
    end
    return s
end
function beta(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[d+i]*t^i
    end
    return s
end
function gamma(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[2*d+i]*t^i
    end
    return s
end
function c1(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[3*d+i]*t^i
    end
    return s
end
function c2(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[4*d+i]*t^i
    end
    return s
end

#mod out global sclaing here
function c3(t,Samp)
    s=0
    for i in (1:d-1)
        s += Samp[5*d+i]*t^i
    end
    s += t^d
    return s
end

function Xv(Samp)
    return [Samp[6*d]; Samp[6*d+1]; Samp[6*d+2]; 1]
end


function P(t,Samp)
    return cayleyNotNorm(alpha(t,Samp), beta(t,Samp), gamma(t,Samp))*[1 0 0 -c1(t,Samp); 0 1 0 -c2(t,Samp); 0 0 1 -c3(t,Samp)]
end

function sl(t)
    return [1 0 -t]
end



function skewmatrix(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end



#############
#Checking the degree of the problems with p=1 using the skew matrix formulation of the problem in [27]



randomSamp=vcat(randn(ComplexF64,2,3*d+1),ones(ComplexF64,1,3*d+1))

eqs=[]
for i in (1:3*d+1)
    push!(eqs, skewmatrix(randomSamp[1:3,i])*cayleyNotNorm(alpha(randomSamp[2,i],x), beta(randomSamp[2,i],x), gamma(randomSamp[2,i],x))*[x[6*d]-c1(randomSamp[2,i],x); x[6*d+1]-c2(randomSamp[2,i],x); x[6*d+2]-c3(randomSamp[2,i],x)])
end

eqns=vcat(eqs...)
F=System(eqns)
solve(F, start_system = :polyhedral)





############################
#The following code was not used in the paper, due unexpected results while applying monodromy. 
#It is added here for the curious reader but it is not fully checked for correctness.


#=


funs=Vector{Expression}(undef, 6*d+2)

for i in (1:3*d+1)
    funs[i]=(sl(U[i])*P(U[i],x)*Xv(x))[1,1]
    funs[3*d+1+i]=U[3*d+1+i]*(P(U[i],x)*Xv(x))[3,1]-(P(U[i],x)*Xv(x))[2,1]
end
F = System(funs; variables = V, parameters = U)

V0 = randn(ComplexF64, 6*d+2)
U0=zeros(ComplexF64, 6*d+2)


f=(sl(t0)*P(t0,x)*Xv(x))[1,1]
F0=System([f], variables = [t0], parameters = V)

S0 = solve(F0, start_system = :total_degree, target_parameters = V0)
sol=solutions(S0)
for i in 1:3*d+1
    U0[i]=(P(sol[i][1],V0)*Xv(V0))[1]/(P(sol[i][1],V0)*Xv(V0))[3]
    #U0[i]=sol[i][1]
    U0[3*d+1+i]=(P(sol[i][1],V0)*Xv(V0))[2]/(P(sol[i][1],V0)*Xv(V0))[3]
end


SM=monodromy_solve(F, V0, U0)




solve(F, start_system = :polyhedral, target_parameters = U0)

#case p=2 d=2 delta=0
#=
@var x[1:22]
@var t0
V=x[1:11]
U=x[12:22]


function P(t)
    return [1 0 0 -x[1]*t-x[2]*t^2; 0 1 0 -x[3]*t-x[4]*t^2; 0 0 1 -x[5]*t-t^2]
end

function Peval(t,Samp)
    return [1 0 0 -Samp[1]*t-Samp[2]*t^2; 0 1 0 -Samp[3]*t-Samp[4]*t^2; 0 0 1 -Samp[5]*t-t^2]
end

X= [x[6]; x[7]; x[8]; 1]
Y= [x[9]; x[10]; x[11]; 1]
function sl(t)
    return [1 0 -t]
end

ratfuns=[(sl(x[12])*P(x[12])*X)[1,1], (sl(x[14])*P(x[14])*X)[1,1], (sl(x[16])*P(x[16])*X)[1,1], (sl(x[18])*P(x[18])*Y)[1,1], (sl(x[20])*P(x[20])*Y)[1,1],
    x[13]*(P(x[12])*X)[3,1]-(P(x[12])*X)[2,1],x[15]*(P(x[14])*X)[3,1]-(P(x[14])*X)[2,1],x[17]*(P(x[16])*X)[3,1]-(P(x[16])*X)[2,1],
    x[19]*(P(x[18])*Y)[3,1]-(P(x[18])*Y)[2,1],x[21]*(P(x[20])*Y)[3,1]-(P(x[20])*Y)[2,1],
    x[22]*(P(x[12]+x[14]+x[16]-x[18]-x[20])*Y)[3,1]-(P(x[12]+x[14]+x[16]-x[18]-x[20])*Y)[2,1]]

F = System(ratfuns; variables = V, parameters = U)
V0 = randn(ComplexF64, 11)

U0=zeros(ComplexF64, 11)
f=(sl(t0)*P(t0)*X)[1,1]
F0=System([f], variables = [t0], parameters = V)
S0 = solve(F0, start_system = :total_degree, target_parameters = V0)
solsX = solutions(S0)

X0=[V0[6]; V0[7]; V0[8]; 1]
Y0=[V0[9]; V0[10]; V0[11]; 1]
U0[1]=solsX[1][1]
U0[3]=solsX[2][1]
U0[5]=solsX[3][1]
for i in 1:3
    U0[2*i]=(Peval(solsX[i][1],V0)*X0)[2]/(Peval(solsX[i][1],V0)*X0)[3]
end
ff=(sl(t0)*P(t0)*Y)[1,1]
F1=System([ff], variables = [t0], parameters = V)
S1 = solve(F1, start_system = :total_degree, target_parameters = V0)
solsY = solutions(S1)
U0[7]=solsY[1][1]
U0[9]=solsY[2][1]
U0[8]=(Peval(solsY[1][1],V0)*Y0)[2]/(Peval(solsY[1][1],V0)*Y0)[3]
U0[10]=(Peval(solsY[2][1],V0)*Y0)[2]/(Peval(solsY[2][1],V0)*Y0)[3]
U0[11]=(Peval(solsY[3][1],V0)*Y0)[2]/(Peval(solsY[3][1],V0)*Y0)[3]
SM=monodromy_solve(F, V0, U0)
=#

=#
