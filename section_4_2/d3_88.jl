using HomotopyContinuation

@var vx,vy,vz, wx,wy,wz, v3x,v3y,v3z, a1,b1,c1,d1, a2,b2,c2
d2 = 1
@var x1a,y1a,x1b,y1b,x1c,y1c,x1d,y1d,x1e,y1e,x1f,y1f,x1g,y1g,x1h,y1h,x2a,y2a,x2b,y2b,x2c,y2c,x2d,y2d,x2e,y2e,x2f,y2f,x2g,y2g,x2h,y2h

variables = [vx,vy,vz, wx,wy,wz, v3x,v3y,v3z, a1,b1,c1,d1, a2,b2,c2]
parameters = [x1a,y1a,x1b,y1b,x1c,y1c,x1d,y1d,x1e,y1e,x1f,y1f,x1g,y1g,x1h,y1h,x2a,y2a,x2b,y2b,x2c,y2c,x2d,y2d,x2e,y2e,x2f,y2f,x2g,y2g,x2h,y2h]

rotRS(a) = [
        1+a[1]*a[1]-a[2]*a[2]-a[3]*a[3]  2*a[1]*a[2]-2*a[3]  2*a[2]+2*a[1]*a[3]
        2*a[3]+2*a[1]*a[2]  1-a[1]*a[1]+a[2]*a[2]-a[3]*a[3]  2*a[2]*a[3]-2*a[1]
        2*a[1]*a[3]-2*a[2]  2*a[1]+2*a[2]*a[3]  1-a[1]*a[1]-a[2]*a[2]+a[3]*a[3]
]

skewMatrix(u) = [
	0 -u[3] u[2]
	u[3] 0 -u[1]
	-u[2] u[1] 0
]

function skewmatrix(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end


L1 = [c1;0;d1]
L2 = [c2;0;d2]
V = [vx;vy;vz]
W = [wx;wy;wz]
V3 = [v3x;v3y;v3z]

Ld1 = [a1;1;b1]
Ld1x = [0 -b1 1; b1 0 -a1; -1 a1 0]
Ld2 = [a2;1;b2]
Ld2x = [0 -b2 1; b2 0 -a2; -1 a2 0]


u1a = [x1a y1a 1]
u1b = [x1b y1b 1]
u1c = [x1c y1c 1]
u1d = [x1d y1d 1]
u1e = [x1e y1e 1]
u1f = [x1f y1f 1]
u1g = [x1g y1g 1]
u1h = [x1h y1h 1]

u2a = [x2a y2a 1]
u2b = [x2b y2b 1]
u2c = [x2c y2c 1]
u2d = [x2d y2d 1]
u2e = [x2e y2e 1]
u2f = [x2f y2f 1]
u2g = [x2g y2g 1]
u2h = [x2h y2h 1]

#u1ax = skewmatrix(u1a)
#u1bx = skewmatrix(u1b)
#u1cx = skewmatrix(u1c)
#u1dx = skewmatrix(u1d)
#u1ex = skewmatrix(u1e)
#u1fx = skewmatrix(u1f)
#u1gx = skewmatrix(u1g)
#u1hx = skewmatrix(u1h)

#u2ax = skewmatrix(u2a)
#u2bx = skewmatrix(u2b)
#u2cx = skewmatrix(u2c)
#u2dx = skewmatrix(u2d)
#u2ex = skewmatrix(u2e)
#u2fx = skewmatrix(u2f)
#u2gx = skewmatrix(u2g)
#u2hx = skewmatrix(u2h)

#R1a = rotRS(x1a*A)
#R1b = rotRS(x1b*A)
#R2a = rotRS(x2a*A)
#R2b = rotRS(x2b*A)
#R3a = rotRS(x3a*A)
#R3b = rotRS(x3b*A)
#R4a = rotRS(x4a*A)
#R4b = rotRS(x4b*A)
#R5a = rotRS(x5a*A)
#R5b = rotRS(x5b*A)

eqn1a = u1a * Ld1x * (L1 - y1a*V - y1a^2*W - y1a^3*V3)
eqn1b = u1b * Ld1x * (L1 - y1b*V - y1b^2*W - y1b^3*V3)
eqn1c = u1c * Ld1x * (L1 - y1c*V - y1c^2*W - y1c^3*V3)
eqn1d = u1d * Ld1x * (L1 - y1d*V - y1d^2*W - y1d^3*V3)
eqn1e = u1e * Ld1x * (L1 - y1e*V - y1e^2*W - y1e^3*V3)
eqn1f = u1f * Ld1x * (L1 - y1f*V - y1f^2*W - y1f^3*V3)
eqn1g = u1g * Ld1x * (L1 - y1g*V - y1g^2*W - y1g^3*V3)
eqn1h = u1h * Ld1x * (L1 - y1h*V - y1h^2*W - y1h^3*V3)

eqn2a = u2a * Ld2x * (L2 - y2a*V - y2a^2*W - y2a^3*V3)
eqn2b = u2b * Ld2x * (L2 - y2b*V - y2b^2*W - y2b^3*V3)
eqn2c = u2c * Ld2x * (L2 - y2c*V - y2c^2*W - y2c^3*V3)
eqn2d = u2d * Ld2x * (L2 - y2d*V - y2d^2*W - y2d^3*V3)
eqn2e = u2e * Ld2x * (L2 - y2e*V - y2e^2*W - y2e^3*V3)
eqn2f = u2f * Ld2x * (L2 - y2f*V - y2f^2*W - y2f^3*V3)
eqn2g = u2g * Ld2x * (L2 - y2g*V - y2g^2*W - y2g^3*V3)
eqn2h = u2h * Ld2x * (L2 - y2h*V - y2h^2*W - y2h^3*V3)


eqns = [eqn1a;eqn1b;eqn1c;eqn1d;eqn1e;eqn1f;eqn1g;eqn1h;eqn2a;eqn2b;eqn2c;eqn2d;eqn2e;eqn2f;eqn2g;eqn2h]
F = System(eqns; variables = variables, parameters = parameters)

#p = rand(length(parameters))
p = randn(ComplexF64,length(parameters))
println(p)
sol = HomotopyContinuation.solve(F;target_parameters=p, show_progress = true)
#sol = monodromy_solve(F;show_progress = true)
print(sol)
