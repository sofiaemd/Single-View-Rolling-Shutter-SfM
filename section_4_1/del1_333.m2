restart
needsPackage "MonodromySolver"
debug MonodromySolver
needsPackage "Posets"


FF = ZZ/32003
R = FF[ax,ay,az, q1x,q1y, q2x,q2y, q3x,q3y]

x1a = random(FF)
y1a = random(FF)
x1b = random(FF)
y1b = random(FF)
x1c = random(FF)
y1c = random(FF)
x2a = random(FF)
y2a = random(FF)
x2b = random(FF)
y2b = random(FF)
x2c = random(FF)
y2c = random(FF)
x3a = random(FF)
y3a = random(FF)
x3b = random(FF)
y3b = random(FF)
x3c = random(FF)
y3c = random(FF)

-*
THE FOLLOWING CODE IS TAKEN FROM PLMP (https://github.com/timduff35/PLMP/blob/master/CODE/common.m2)
-- produces gates for "small" determinants"
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M_{1,2}^{1,2})-M_(0,1)*det2(M_{0,2}^{1,2})+M_(0,2)*det2(M_{0,1}^{1,2})

END OF THE PLMP CODE
*-

rotRS = a -> matrix{
        {1+a_(0,0)*a_(0,0)-a_(1,0)*a_(1,0)-a_(2,0)*a_(2,0), 2*a_(0,0)*a_(1,0)-2*a_(2,0), 2*a_(1,0)+2*a_(0,0)*a_(2,0)},
        {2*a_(2,0)+2*a_(0,0)*a_(1,0), 1-a_(0,0)*a_(0,0)+a_(1,0)*a_(1,0)-a_(2,0)*a_(2,0), 2*a_(1,0)*a_(2,0)-2*a_(0,0)},
        {2*a_(0,0)*a_(2,0)-2*a_(1,0), 2*a_(0,0)+2*a_(1,0)*a_(2,0), 1-a_(0,0)*a_(0,0)-a_(1,0)*a_(1,0)+a_(2,0)*a_(2,0)}
    }


A = matrix{{ax},{ay},{az}}

---

q1 = matrix{{q1x},{q1y},{1}}

R1a = rotRS(y1a*A)
r1ax = matrix{{0,y1a,1},{-y1a,0,0},{-1,0,0}}
u1ax = matrix{{0,-1,y1a},{1,0,-x1a},{-y1a,x1a,0}}

eqn1a = u1ax * r1ax * R1a * q1

R1b = rotRS(y1b*A)
r1bx = matrix{{0,y1b,1},{-y1b,0,0},{-1,0,0}}
u1bx = matrix{{0,-1,y1b},{1,0,-x1b},{-y1b,x1b,0}}

eqn1b = u1bx * r1bx * R1b * q1

R1c = rotRS(y1c*A)
r1cx = matrix{{0,y1c,1},{-y1c,0,0},{-1,0,0}}
u1cx = matrix{{0,-1,y1c},{1,0,-x1c},{-y1c,x1c,0}}

eqn1c = u1cx * r1cx * R1c * q1



q2 = matrix{{q2x},{q2y},{1}}

R2a = rotRS(y2a*A)
r2ax = matrix{{0,y2a,1},{-y2a,0,0},{-1,0,0}}
u2ax = matrix{{0,-1,y2a},{1,0,-x2a},{-y2a,x2a,0}}

eqn2a = u2ax * r2ax * R2a * q2

R2b = rotRS(y2b*A)
r2bx = matrix{{0,y2b,1},{-y2b,0,0},{-1,0,0}}
u2bx = matrix{{0,-1,y2b},{1,0,-x2b},{-y2b,x2b,0}}

eqn2b = u2bx * r2bx * R2b * q2

R2c = rotRS(y2c*A)
r2cx = matrix{{0,y2c,1},{-y2c,0,0},{-1,0,0}}
u2cx = matrix{{0,-1,y2c},{1,0,-x2c},{-y2c,x2c,0}}

eqn2c = u2cx * r2cx * R2c * q2


q3 = matrix{{q3x},{q3y},{1}}

R3a = rotRS(y3a*A)
r3ax = matrix{{0,y3a,1},{-y3a,0,0},{-1,0,0}}
u3ax = matrix{{0,-1,y3a},{1,0,-x3a},{-y3a,x3a,0}}

eqn3a = u3ax * r3ax * R3a * q3

R3b = rotRS(y3b*A)
r3bx = matrix{{0,y3b,1},{-y3b,0,0},{-1,0,0}}
u3bx = matrix{{0,-1,y3b},{1,0,-x3b},{-y3b,x3b,0}}

eqn3b = u3bx * r3bx * R3b * q3

R3c = rotRS(y3c*A)
r3cx = matrix{{0,y3c,1},{-y3c,0,0},{-1,0,0}}
u3cx = matrix{{0,-1,y3c},{1,0,-x3c},{-y3c,x3c,0}}

eqn3c = u3cx * r3cx * R3c * q3



eqns = eqn1a^{1} || eqn1b^{1} || eqn1c^{1} || eqn2a^{1} || eqn2b^{1} || eqn2c^{1} || eqn3a^{1} || eqn3b^{1} || eqn3c^{1}

Ispec = ideal eqns
elapsedTime G := groebnerBasis(Ispec, Strategy => "F4");
ltG := ideal leadTerm G; (dim ltG, degree ltG);
dimdeg := (dim ltG, degree ltG);
<< "(dim, degree) " << dimdeg << endl;


