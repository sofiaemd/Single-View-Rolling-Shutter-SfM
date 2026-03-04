restart

d = 2
l = 6
R = QQ[c_(0,1) .. c_(2,d),q_0, q_1,q_2,D_0,D_1,D_2][x]

cP = (u,w) -> (
u = flatten entries u;
w = flatten entries w;
matrix{{u_1*w_2-u_2*w_1},{u_2*w_0-u_0*w_2},{u_0*w_1-u_1*w_0}}
)

randomPluck = () -> (
y1 = apply(3, i -> random(-20,20));
y2 = apply(3, i -> random(-20,20));
pl = flatten(entries(cP(matrix({y1}),matrix({y2})))) | flatten entries (vector(y1) - vector(y2));
apply(pl, v -> v/pl_2)
)

evalRandomLine = () -> (
pl = randomPluck();
{q_0 => pl_0, q_1 => pl_1, q_2 => pl_2, D_0 => pl_3, D_1 => pl_4, D_2 => pl_5}
)

evalRandomLineI = (i) -> (
pl = randomPluck();
{q_(i,0) => pl_0, q_(i,1) => pl_1, D_(i,0) => pl_3, D_(i,1) => pl_4}
)
    
evalRandomCamera = () -> (
randomCoeffs = apply(3*d, i->random(-20,20));
normCoeffs = apply(randomCoeffs, c->c/randomCoeffs_(3*d-1));
var = toList(c_(0,1) .. c_(2,d));
apply(3*d, i -> var_i => normCoeffs_i)
)


C = matrix({apply(3, j-> sum(apply (d, i->c_(j,i+1)*x^(i+1))))})

Q = transpose matrix({{q_0,q_1,q_2}})
Delta = transpose matrix({{D_0,D_1,D_2}})
r = transpose matrix({{1,0,-x}})

lineCurve = cP(r, Q + cP(Delta,C))


-- testing the l lines
S1 = QQ[c_(0,1) .. c_(2,d),q_0, q_1,q_2,D_0,D_1,D_2,
    q_(0,0) .. q_(l-1,0), q_(0,1) .. q_(l-1,1),
    D_(0,0) .. D_(l-1,0), D_(0,1) .. D_(l-1,1)]
coeffs = apply(d+1, i -> sub((entries last(coefficients((entries lineCurve)_0_0)))_i_0,S1))|
apply(d+2, i -> sub((entries last(coefficients((entries lineCurve)_1_0)))_i_0,S1))

interestingCoeffs = coeffs_({0,1,2,3,5})

evalCoeffs = apply(l, i -> {q_0 => q_(i,0), q_1 => q_(i,1), q_2 => 1, 
        D_0 => D_(i,0), D_1 => D_(i,1), D_2 => -q_(i,0)*D_(i,0)-q_(i,1)*D_(i,1)})

evalCam = evalRandomCamera()
evalLine = apply(l, i -> evalRandomLine())


eqs1 = flatten apply(l, i -> apply(interestingCoeffs, p -> sub(p, evalCoeffs_i|{c_(2,d) => 1}) - sub(p, evalCam|evalLine_i)))
eqs0 = flatten apply(l, i -> apply(interestingCoeffs, p -> sub(p, evalCoeffs_i|{c_(2,d) => 0}) - sub(p, evalCam|evalLine_i)))

S2 = QQ[c_(0,1) .. c_(2,d-1),c_(0,d),c_(1,d),
    q_(0,0) .. q_(l-1,0), q_(0,1) .. q_(l-1,1),
    D_(0,0) .. D_(l-1,0), D_(0,1) .. D_(l-1,1)]

I = ideal(apply(eqs1, p -> sub(p,S2)))
dim(I)
degree(I)

J = ideal(apply(eqs0, p -> sub(p,S2)))
dim(J)
degree(J)

-- jacobian

linesToCurves = flatten apply(l, i -> apply(interestingCoeffs, p -> sub(sub(p,evalCoeffs_i|{c_(2,d) => 1}),S2)))
j = jacobian(matrix({linesToCurves_(toList{0..28})}))
det(sub(j,flatten(evalCam_(toList{0..4})|apply(l, i -> evalRandomLineI(i)))))





            
            
            
            
            
            
