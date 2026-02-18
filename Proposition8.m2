restart

d = 2
delta = 2
l = 3
S = QQ[c_(0,1) .. c_(2,d),a_(0,1) .. a_(2,delta),
    q_1,q_2,D_1,D_2,
    q_(0,1)..q_(l-1,1),q_(0,2)..q_(l-1,2),D_(0,1)..D_(l-1,1),D_(0,2)..D_(l-1,2)]
R = S[x]

cP = (u,w) -> (
u = flatten entries u;
w = flatten entries w;
matrix{{u_1*w_2-u_2*w_1},{u_2*w_0-u_0*w_2},{u_0*w_1-u_1*w_0}}
)

cayleyTransform = () -> (
A = apply(3, j-> sum(apply (delta, i->a_(j,i+1)*x^(i+1))));
matrix({{1+A_0^2-A_1^2-A_2^2,2*(A_0*A_1-A_2),2*(A_0*A_2+A_1)},
        {2*(A_0*A_1+A_2),1-A_0^2+A_1^2-A_2^2, 2*(A_1*A_2-A_0)},
        {2*(A_0*A_2-A_1),2*(A_1*A_2+A_0),1-A_0^2-A_1^2+A_2^2}})
)

C = matrix({apply(3, j-> sum(apply (d, i->c_(j,i+1)*x^(i+1))))})

-- PASO1: calcular los coeffs que nos interesan

Q = transpose matrix({{q_1,q_2,1}})
Delta = transpose matrix({{D_1,D_2,-q_1*D_1-q_2*D_2}})
r = transpose matrix({{1,0,-x}})

lineCurve = cP(r, cayleyTransform()*(Q + cP(Delta,C)))

coeffs = (apply(2*delta+d+1, i -> sub((entries last(coefficients((entries lineCurve)_0_0)))_i_0,S)))|
(apply(2*delta+d+2, i -> sub((entries last(coefficients((entries lineCurve)_1_0)))_i_0,S)))
coeffs = apply(coeffs, p -> sub(p,S))

--PASO2: reemplazar por coeffs random y generar ideal

evalCam = apply(splice toList(c_(0,1) .. c_(2,d-1), c_(0,d),c_(1,d),a_(0,1) .. a_(2,delta)), 
    v -> v => random(-20,20))
evalCam1 = evalCam|{c_(2,d) => 1}
evalCam0 = evalCam|{c_(2,d) => 0}
-- evalLines = apply(l, i -> {q_1 => random(-20,20), q_2 => random(-20,20), 
--        D_1 => random(-20,20), D_2 => random(-20,20)})
evalLines = {{q_1 => 0, q_2 => 0, D_1 => 1, D_2 => 0},
    {q_1 => 0, q_2 => 0, D_1 => 0, D_2 => 1},
    {q_1 => 1, q_2 => 0, D_1 => 1, D_2 => 0}}
evalCoeffs = apply(l, i -> {q_1 => q_(i,1), q_2 => q_(i,2), 
        D_1 => D_(i,1), D_2 => D_(i,2)})

S1 = QQ[c_(0,1) .. c_(2,d-1), c_(0,d),c_(1,d),a_(0,1) .. a_(2,delta),
    q_1,q_2,D_1,D_2,
    q_(0,1)..q_(l-1,1),q_(0,2)..q_(l-1,2),D_(0,1)..D_(l-1,1),D_(0,2)..D_(l-1,2)]
    
eqs1 = flatten apply(l, i -> apply(coeffs, p -> sub(sub(p, evalCoeffs_i) - sub(p, evalCam1|evalLines_i),S1)))
eqs0 = flatten apply(l, i -> apply(coeffs, p -> sub(sub(p, evalCoeffs_i) - sub(p, evalCam0|evalLines_i),S1)))

I=ideal(eqs1)
dim(I)
degree(I)

J = ideal(eqs0)
dim(J)
degree(J)

-- Jacobian

jacMatrix = jacobian(matrix({eqs1}))
jacDet = det(submatrix(jacMatrix, ,splice toList{1..2*delta,2*delta+2..3*delta+3}))

sub(jacDet, eval_(toList(splice 0..3*delta+1)))

