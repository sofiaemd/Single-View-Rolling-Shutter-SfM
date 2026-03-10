restart

-- Defining the variables and polynomial ring (here you can choose d = 1, 2 or 3)
d = 1
delta = 1
S = QQ[c_(0,1) .. c_(0,d),c_(2,1)..c_(2,d-1),a_(0,1) .. a_(2,delta),
    q_1,q_2,D_1,D_2]
R = S[x]
c_(2,d)=1
D_3=-q_1*D_1-q_2*D_2

for i from 1 to d do(c_(1,i)=0)

-- Function that computes cross product
cP = (u,w) -> (
u = flatten entries u;
w = flatten entries w;
matrix{{u_1*w_2-u_2*w_1},{u_2*w_0-u_0*w_2},{u_0*w_1-u_1*w_0}}
)

-- Functions that returns the rotation matrix
cayleyTransform = () -> (
A = apply(3, j-> sum(apply (delta, i->a_(j,i+1)*x^(i+1))));
matrix({{1+A_0^2-A_1^2-A_2^2,2*(A_0*A_1-A_2),2*(A_0*A_2+A_1)},
        {2*(A_0*A_1+A_2),1-A_0^2+A_1^2-A_2^2, 2*(A_1*A_2-A_0)},
        {2*(A_0*A_2-A_1),2*(A_1*A_2+A_0),1-A_0^2-A_1^2+A_2^2}})
)

-- Creating the polynomials in the camera movement
C = matrix({apply(3, j-> sum(apply (d, i->c_(j,i+1)*x^(i+1))))})


-- Computing polynomials defining the picture-taking map
Q = transpose matrix({{q_1,q_2,1}})
Delta = transpose matrix({{D_1,D_2,D_3}})
r = transpose matrix({{1,0,-x}})

lineCurve = cP(r, cayleyTransform()*(Q + cP(Delta,C)))

coeffs = (apply(2*delta+d+1, i -> sub((entries last(coefficients((entries lineCurve)_1_0)))_i_0,S)))|
(apply(2*delta+d+1, i -> sub((entries last(coefficients((entries lineCurve)_0_0)))_i_0,S)))

-- Computing the Jacobian for the picture-taking map
J=jacobian(matrix{coeffs})

-- Choosing a point to evaluate
if (d>1) then evalParams = apply(d, i -> c_(0,i+1) => 0)|{c_(2,1)=>1/2}|apply(d-2, i -> c_(2,i+2) => 0)|
{q_1=>0,q_2=>1,D_1=>0,D_2=>1,a_(0,1)=>0,a_(1,1)=>1/2,a_(2,1)=>0}
if d==1 then evalParams = {c_(0,1)=>0,q_1=>0,q_2=>1,D_1=>0,D_2=>1,a_(0,1)=>0,a_(1,1)=>1/2,a_(2,1)=>0}

-- Evaluating the jacobian and computing the determinant
det(sub(J,evalParams))