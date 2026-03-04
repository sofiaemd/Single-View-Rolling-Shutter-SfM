using LinearAlgebra

x1,y1 = 1,2
x2,y2 = 3,5
x3,y3 = 7,11
x4,y4 = 13,17
x5,y5 = 19,23
y6  = 31
z=x1+x2+x3-x4-x5

M = [
 -x1 0 x1^2 1 0 -x1 0 0 0;
 0 -x1 -y1*x1 0 -1 y1 0 0 0;
 -x2 0 x2^2 1 0 -x2 0 0 0;
 0 -x2 -y2*x2 0 -1 y2 0 0 0;
 -x3 0 x3^2 0 0 0 1 0 -x3;
 0 -x3 -y3*x3 0 0 0 0 -1 y3;
 -x4 0 x4^2 0 0 0 1 0 -x4;
 0 -x4 -y4*x4 0 0 0 0 -1 y4
]

N = [
 -x1  -x1^2   0     0    x1^2   x1^3   1   0  -x1   0   0   0;
  0     0    -x1  -x1^2  -y1*x1 -y1*x1^2 0  -1   y1   0   0   0;
 -x2  -x2^2   0     0    x2^2   x2^3   1   0  -x2   0   0   0;
  0     0    -x2  -x2^2  -y2*x2 -y2*x2^2 0  -1   y2   0   0   0;
 -x3  -x3^2   0     0    x3^2   x3^3   1   0  -x3   0   0   0;
  0     0    -x3  -x3^2  -y3*x3 -y3*x3^2 0  -1   y3   0   0   0;
 -x4  -x4^2   0     0    x4^2   x4^3   0   0    0    1   0  -x4;
  0     0    -x4  -x4^2  -y4*x4 -y4*x4^2 0   0    0    0  -1   y4;
 -x5  -x5^2   0     0    x5^2   x5^3   0   0    0    1   0  -x5;
  0     0    -x5  -x5^2  -y5*x5 -y5*x5^2 0   0    0    0  -1   y5;
  0     0    -z   -z^2   -y6*z  -y6*z^2  0   0    0    0  -1   y6
]
d=3
p=3*d-1

M2=zeros(ComplexF64, p,p)
for i in (1:p)
    rs=randn(ComplexF64, 3)
    for j in (1:d)
        M2[i,j]=(rs[3]-rs[2])*rs[1]^(j-1)
    end
    for j in (1:d)
        M2[i,j+d]=rs[1]^(j)
    end
    for j in (1:d-1)
        M2[i,j+2*d]=-rs[3]*rs[1]^(j)
    end
end
