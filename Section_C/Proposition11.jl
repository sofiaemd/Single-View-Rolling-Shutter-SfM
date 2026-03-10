using HomotopyContinuation
using LinearAlgebra

#Choose here the parameters you want to test.
d=1
delta=3
p=3

@var v,z

C1 = randn(ComplexF64, d)
C2 = randn(ComplexF64, d)
C3 = randn(ComplexF64, d)
C4 = randn(ComplexF64, delta)
C5 = randn(ComplexF64, delta)
C6 = randn(ComplexF64, delta)

#Define the polynommials from the coefficients
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

#Cayley matrix without normalization factor
function cayleyNotNorm(a,b,c)
    fr=1+a^2+b^2+c^2
    M=zeros(typeof(a),3,3)
    M[1,1]=1+a^2-b^2-c^2
    M[1,2]= 2*(a*b-c)
    M[1,3]= 2*(a*c+b)
    M[2,1]=2*(a*b+c)
    M[2,2]= 1-a^2+b^2-c^2
    M[2,3]= 2*(b*c-a)
    M[3,1]=2*(a*c-b)
    M[3,2]= 2*(b*c+a)
    M[3,3]=1-a^2-b^2+c^2
    return M
end 

#Camera matrix

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

#Define scanline

if d==0
    sl=zeros(typeof(r1(v)),1,3)
else
    sl=zeros(typeof(c1(v)),1,3)
end
sl[1]=1
sl[3]= -v
sleq= sl*PmatrixVar(r1,r2,r3,c1,c2,c3)

#Fill matrix A with the images of the random points under a fixed camera 
#The image will live in an affine hyperplane, similalry to Lemma1 if d>1
#but we could vary the camera to avoid this .

A=zeros(ComplexF64, 2*p*(1+d+2*delta), 2*p*(1+d+2*delta))

for i in (1:2*p*(1+d+2*delta))  
    X = randn(ComplexF64, 4,p)
    k=1
    while k <= p
        f=X[1,k]*sleq[1] + X[2,k]*sleq[2] + X[3,k]*sleq[3] + X[4,k]*sleq[4]
        F = System([f, z*(1+r1(v)^2 + r2(v)^2 + r3(v)^2) -1])
        result = HomotopyContinuation.solve(F)
        sols=solutions(result)
        if length(sols)<(1+d+2*delta) 
            #catching numerical issues when solving the equation as random coeefs are small
            #println("Not enough solutions found!")
        else
            for j in (1:(1+d+2*delta))
                P=PmatrixVal(r1,r2,r3,c1,c2,c3,sols[j][1])
                im= P*X[1:4,k]
                A[i,1+2*(j-1)+(k-1)*2*(1+d+2*delta)]=im[1]/im[3]
                A[i,2+2*(j-1)+(k-1)*2*(1+d+2*delta)]=im[2]/im[3]
            end
            k += 1
        end
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
println("Numerical rank of the constructed matrix A is ", numranktest(A))
println("Maximal rank is ", 2*p*(1+d+2*delta))
if d>1
    println("Expected rank drop due to Lemma 1 is ", p-1)
    if numranktest(A)+p-1== 2*p*(1+d+2*delta)
        println("Proposition 11 is numerically verified in this case.")
    else
        println("Proposition 11 is not numerically verified, potentially choose bigger coefficients.")
    end
else
    if numranktest(A)== 2*p*(1+d+2*delta)
        println("Proposition 11 is numerically verified in this case.")
    else
        println("Proposition 11 is not numerically verified, potentially choose bigger coefficients.")
    end
end



