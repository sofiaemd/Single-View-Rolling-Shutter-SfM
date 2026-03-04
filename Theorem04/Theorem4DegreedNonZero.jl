using HomotopyContinuation

d=1
delta=1
p=3*(d+delta)-1
@var x[1:4*p]
#@var t0



function cayleyNotNorm(a,b,c)
    return [1+a^2-b^2-c^2 2*(a*b-c) 2*(a*c+b); 2*(a*b+c) 1-a^2+b^2-c^2 2*(b*c-a); 2*(a*c-b) 2*(b*c+a) 1-a^2-b^2+c^2]
end

function alpha(t,Samp)
    s=0
    for i in (1:delta)
        s += Samp[i]*t^i
    end
    return s
end
function beta(t,Samp)
    s=0
    for i in (1:delta)
        s += Samp[delta+i]*t^i
    end
    return s
end
function gamma(t,Samp)
    s=0
    for i in (1:delta)
        s += Samp[2*delta+i]*t^i
    end
    return s
end
function c1(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[3*delta+i]*t^i
    end
    return s
end
function c2(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[3*delta+d+i]*t^i
    end
    return s
end
function c3(t,Samp)
    s=0
    for i in (1:d)
        s += Samp[3*delta+2*d+i]*t^i
    end
    return s
end

Xs=vcat(zeros(Expression,3,p),ones(Expression,1,p))
for i in (1:p-1)
    Xs[1,i]=x[p+1+3*(i-1)+1]
    Xs[2,i]=x[p+1+3*(i-1)+2]
    Xs[3,i]=x[p+1+3*(i-1)+3]
end
Xs[1,p]=x[4*p-1]
Xs[2,p]=x[4*p]
Xs[3,p]=1

function P(t,Samp)
    return cayleyNotNorm(alpha(t,Samp), beta(t,Samp), gamma(t,Samp))*[1 0 0 -c1(t,Samp); 0 1 0 -c2(t,Samp); 0 0 1 -c3(t,Samp)]
end

function sls(t)
    return [1 0 -t]
end

function skewmatrix(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function cayleyNotNorm(a,b,c)
    return [1+a^2-b^2-c^2 2*(a*b-c) 2*(a*c+b); 2*(a*b+c) 1-a^2+b^2-c^2 2*(b*c-a); 2*(a*c-b) 2*(b*c+a) 1-a^2-b^2+c^2]
end

#randomSamp=vcat(randn(ComplexF64,2,2*p),ones(ComplexF64,1,2*p))


##########################
#Skew matrix formulation



eqs=[]
for i in (1:p)
    local randomSamp=vcat(randn(ComplexF64,2,2),ones(ComplexF64,1,2))
    push!(eqs, skewmatrix(randomSamp[1:3,1])*cayleyNotNorm(alpha(randomSamp[2,1],x), beta(randomSamp[2,1],x), gamma(randomSamp[2,1],x))*[Xs[1,i]-c1(randomSamp[2,1],x); Xs[2,i]-c2(randomSamp[2,1],x); Xs[3,i]-c3(randomSamp[2,1],x)])
    push!(eqs, skewmatrix(randomSamp[1:3,2])*cayleyNotNorm(alpha(randomSamp[2,2],x), beta(randomSamp[2,2],x), gamma(randomSamp[2,2],x))*[Xs[1,i]-c1(randomSamp[2,2],x); Xs[2,i]-c2(randomSamp[2,2],x); Xs[3,i]-c3(randomSamp[2,2],x)])

end

eqns=vcat(eqs...)
F=System(eqns)
S0=solve(F, start_system = :polyhedral)

