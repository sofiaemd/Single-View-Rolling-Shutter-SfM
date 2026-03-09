restart
needsPackage "MonodromySolver"
debug MonodromySolver
needsPackage "Posets"


FF = ZZ/32003
R = FF[ax,ay,az, q1x,q1y]

x1a = random(FF)
y1a = random(FF)
x1b = random(FF)
y1b = random(FF)
x1c = random(FF)
y1c = random(FF)
x1d = random(FF)
y1d = random(FF)
x1e = random(FF)
y1e = random(FF)

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

R1d = rotRS(y1d*A)
r1dx = matrix{{0,y1d,1},{-y1d,0,0},{-1,0,0}}
u1dx = matrix{{0,-1,y1d},{1,0,-x1d},{-y1d,x1d,0}}

eqn1d = u1dx * r1dx * R1d * q1

R1e = rotRS(y1e*A)
r1ex = matrix{{0,y1e,1},{-y1e,0,0},{-1,0,0}}
u1ex = matrix{{0,-1,y1e},{1,0,-x1e},{-y1e,x1e,0}}

eqn1e = u1ex * r1ex * R1e * q1


eqns = eqn1a^{1} || eqn1b^{1} || eqn1c^{1} || eqn1d^{1} || eqn1e^{1}

Ispec = ideal eqns
elapsedTime G := groebnerBasis(Ispec, Strategy => "F4");
ltG := ideal leadTerm G; (dim ltG, degree ltG);
dimdeg := (dim ltG, degree ltG);
<< "(dim, degree) " << dimdeg << endl;


