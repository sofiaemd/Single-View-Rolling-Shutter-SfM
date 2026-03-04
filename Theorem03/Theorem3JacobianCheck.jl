using HomotopyContinuation
using LinearAlgebra


d=5


p=1
@var v,z


delta=d
#We use the special point used in the paper.
X = randn(ComplexF64, 4,p)
X[1]=0
X[2]=6
X[3]=7.5
X[4]=1

C2=zeros(ComplexF64, d)
C3=zeros(ComplexF64, d)
C1=zeros(ComplexF64, d)
C4=zeros(ComplexF64, delta)
C4[delta]=1
C5=zeros(ComplexF64, delta)
C6=zeros(ComplexF64, delta)

C3[d]=1


#The functions here are the same as in the code for Proposition 11.

function c1(t)
    s=0
    for i in (1:d)
        s += C1[i]*t^i
    end
    return s
end
function c2(t)
    s=0
    for i in (1:d)
        s += C2[i]*t^i
    end
    return s
end
function c3(t)
    s=0
    for i in (1:d)
        s += C3[i]*t^i
    end
    return s
end

function r1(t)
    s=0
    for i in (1:delta)
        s += C4[i]*t^i
    end
    return s
end
function r2(t)
    s=0
    for i in (1:delta)
        s += C5[i]*t^i
    end
    return s
end
function r3(t)
    s=0
    for i in (1:delta)
        s += C6[i]*t^i
    end
    return s
end


function PmatrixVar(rr1,rr2,rr3,cc1,cc2,cc3)
    R=cayleyNotNorm(rr1(v),rr2(v),rr3(v))
    if d==0
        P=zeros(typeof(rr1(v)),3,4)
    else
        P=zeros(typeof(cc1(v)),3,4)
    end
    P[1,1]=1
    P[2,2]=1
    P[3,3]=1
    P[1,4]=-cc1(v)
    P[2,4]=-cc2(v)
    P[3,4]=-cc3(v)
    return R*P
end

function PmatrixVal(rr1,rr2,rr3,cc1,cc2,cc3,t)
    R=cayleyNotNorm(ComplexF64(rr1(t)),ComplexF64(rr2(t)),ComplexF64(rr3(t)))
    if d==0
        P=zeros(typeof(rr1(t)),3,4)
    else
        P=zeros(typeof(cc1(t)),3,4)
    end
    P[1,1]=1
    P[2,2]=1
    P[3,3]=1
    P[1,4]=-cc1(t)
    P[2,4]=-cc2(t)
    P[3,4]=-cc3(t)
    return R*P
end


if d==0
    sline=zeros(typeof(r1(v)),1,3)
else
    sline=zeros(typeof(c1(v)),1,3)
end
sline[1]=1
sline[3]= -v
sleq= sline*PmatrixVar(r1,r2,r3,c1,c2,c3)

#The derivatives are computed via the formulas given in the paper.

function delalpha(t,i::Integer)
    M=zeros(typeof(t),3,3)
    M[1,1]=r1(t)
    M[1,2]=r2(t)
    M[1,3]=r3(t)
    M[2,1]=r2(t)
    M[2,2]=-r1(t)
    M[2,3]=-1
    M[3,1]=r3(t)
    M[3,2]=1
    M[3,3]=-r1(t)
    P=zeros(typeof(t),3,4)
    P[1,1]=1
    P[2,2]=1
    P[3,3]=1
    P[1,4]=-c1(t)
    P[2,4]=-c2(t)
    P[3,4]=-c3(t)
    return 2*t^i*M*P*X
end

function delbeta(t,i::Integer)
    M=zeros(typeof(t),3,3)
    M[1,1]=-r2(t)
    M[1,2]=r1(t)
    M[1,3]=1
    M[2,1]=r1(t)
    M[2,2]=r2(t)
    M[2,3]=r3(t)
    M[3,1]=1
    M[3,2]=r3(t)
    M[3,3]=-r2(t)
    P=zeros(typeof(t),3,4)
    P[1,1]=1
    P[2,2]=1
    P[3,3]=1
    P[1,4]=-c1(t)
    P[2,4]=-c2(t)
    P[3,4]=-c3(t)
    return 2*t^i*M*P*X
end

function delgamma(t,i::Integer)
    M=zeros(typeof(t),3,3)
    M[1,1]=-r3(t)
    M[1,2]=-1
    M[1,3]=r1(t)
    M[2,1]=1
    M[2,2]=-r3(t)
    M[2,3]=r2(t)
    M[3,1]=r3(t)
    M[3,2]=r1(t)
    M[3,3]=r3(t)
    P=zeros(typeof(t),3,4)
    P[1,1]=1
    P[2,2]=1
    P[3,3]=1
    P[1,4]=-c1(t)
    P[2,4]=-c2(t)
    P[3,4]=-c3(t)
    return 2*t^i*M*P*X
end

function dela(t,i::Integer)
    M=cayleyNotNorm(r1(t),r2(t),r3(t))
    P=zeros(typeof(t),3,4)
    P[1,4]=-t^i
    return M*P*X
end

function delb(t,i::Integer)
    M=cayleyNotNorm(r1(t),r2(t),r3(t))
    P=zeros(typeof(t),3,4)
    P[2,4]=-t^i
    return M*P*X
end

function delc(t,i::Integer)
    M=cayleyNotNorm(r1(t),r2(t),r3(t))
    P=zeros(typeof(t),3,4)
    P[3,4]=-t^i
    return M*P*X
end

function delx(t,i::Integer)
    P=PmatrixVal(r1,r2,r3,c1,c2,c3,t)
    xx=zeros(typeof(t),4,1)
    xx[i]=1
    return P*xx
end

#We stack them all together to get the curly J matrix which was the Jacobian of P(t)X

function JacP(t)
    J=zeros(typeof(t),3, 3*(d+delta)+2)
    for i in (1:delta)
        J[:,i]=delalpha(t,i)
        J[:,i+delta]=delbeta(t,i)
        J[:,i+2*delta]=delgamma(t,i)
    end
    for i in (1:d)
        J[:,i+3*delta]=dela(t,i)
        J[:,i+3*delta + d]=delb(t,i)
    end
    for i in (1:d-1)
        J[:,i+3*delta + 2*d]=delc(t,i)
    end
    for i in (1:3)
        J[:,i+3*delta + 3*d-1]=delx(t,i)
    end
    return J
end

ff=X[1]*sleq[1] + X[2]*sleq[2] + X[3]*sleq[3] + X[4]*sleq[4]
FF = System([ff, z*(1+r1(v)^2 + r2(v)^2 + r3(v)^2) -1])
result = HomotopyContinuation.solve(FF)
sol=solutions(result)



M= zeros(typeof(sol[1][1]), 2*(1+d+2*delta), 2*(1+d+2*delta))
if length(sol)<(1+d+2*delta)
    println("Not enough solutions found!")
else

    #Now we use the derived formulas to get the Jacobian of the imaging maps which used the implicit functions.    
    for i in 1:(1+d+2*delta)
        J=JacP(sol[i][1])
        slt=zeros(typeof(sol[i][1]), 1,3)
        slt[1]=1
        slt[3]= -sol[i][1]
        e2=zeros(typeof(sol[i][1]), 1,3)
        e2[2]=1
        e3=zeros(typeof(sol[i][1]), 1,3)
        e3[3]=1
        M[:,i]=slt*J
        Px=PmatrixVal(r1,r2,r3,c1,c2,c3,sol[i][1])*X
        M[:,(1+d+2*delta)+i]=e2*J*(Px[3,1])-e3*J*(Px[2,1])
    end
end


function numranktest(M)
    S=svd(M)
    r=0
    for s in S.S
        if s>1e-10
            r+=1
        end
    end
    return r
end

println("The rank of the Jacobian is: ", numranktest(M))
println("The maximal rank is: ", 2*(1+d+2*delta))
if numranktest(M)<2*(1+d+2*delta)
    println("The Jacobian is not full rank, the problem is probably not minimal.")
else
    println("The Jacobian is numerically full rank, the problem is probably minimal.")
end

