restart

-- here you can choose delta
delta = 2
R = QQ[a_(0,1) .. a_(2,delta)][q_0, q_1,q_2][x]

-- Function that computes the cross product
cP = (u,w) -> (
u = flatten entries u;
w = flatten entries w;
matrix{{u_1*w_2-u_2*w_1},{u_2*w_0-u_0*w_2},{u_0*w_1-u_1*w_0}}
)

-- Function that computes the rotation matrix
cayleyTransform = () -> (
b = apply(3, j-> sum(apply (delta, i->a_(j,i+1)*x^(i+1))));
(1/(1+b_0^2+b_1^2+b_2^2))*matrix({{1+b_0^2-b_1^2-b_2^2,2*(b_0*b_1-b_2),2*(b_0*b_2+b_1)},
        {2*(b_0*b_1+b_2),1-b_0^2+b_1^2-b_2^2, 2*(b_1*b_2-b_0)},
        {2*(b_0*b_2-b_1),2*(b_1*b_2+b_0),1-b_0^2-b_1^2+b_2^2}})
)


-- Function that computes the polynomials of the picture-taking map
lineCurveCoeffs = () -> (
R = QQ[a_(0,1) .. a_(2,delta)][q_0, q_1,q_2][x];
L = transpose matrix({{q_0,q_1,q_2}});
r = transpose matrix({{1,0,-x}});
C = cP(r,cayleyTransform()*L);
S = QQ[a_(0,1) .. a_(2,delta), q_0, q_1, q_2];
apply(flatten entries((coefficients(numerator(C_0_1/C_0_2)))#1), c -> sub(c,S)) | 
apply(flatten entries((coefficients(denominator(C_0_1/C_0_2)))#1), c -> sub(c,S))
)

-- Choosing a random camera P to evaluate
eval = apply(toList(a_(0,1) .. a_(2,delta))|{q_0,q_1}, var -> var => random(-20,20))|{q_2 => 1}

-- Computing the ideal that defines the variety of points that have
-- the same image through the picture-taking map as P
I = ideal(apply(lineCurveCoeffs(), p -> p - sub(p, eval)))

-- Verifying the dimension is 0 and degree is 1
dim(I)
degree(I)

-- Performing the jacobian check
T = QQ[a_(0,1) .. a_(2,delta), q_0, q_1]
jacMatrix = jacobian(matrix({apply(lineCurveCoeffs(), p -> sub(p,T))}))
jacDet = det(submatrix(jacMatrix, ,splice toList{1..2*delta,2*delta+2..3*delta+3}))

-- Evaluating jacobian and verifying its determinant is not 0
sub(jacDet, eval_(toList(splice 0..3*delta+1)))
