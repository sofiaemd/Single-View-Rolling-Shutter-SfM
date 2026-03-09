using HomotopyContinuation

@var vx,vy,vz, wx,wy,wz, a1,b1,c1,d1, a2,b2,c2,d2, a3,b3,c3,d3, a4,b4,c4,d4, a5,b5,c5
d5 = 1
@var x1a, y1a, x1b, y1b, x1c, y1c, x1d, y1d, x1e, y1e, x2a, y2a, x2b, y2b, x2c, y2c, x2d, y2d, x2e, y2e, x3a, y3a, x3b, y3b, x3c, y3c, x3d, y3d, x3e, y3e, x4a, y4a, x4b, y4b, x4c, y4c, x4d, y4d, x4e, y4e, x5a, y5a, x5b, y5b, x5c, y5c, x5d, y5d, x5e, y5e

variables = [vx,vy,vz, wx,wy,wz, a1,b1,c1,d1, a2,b2,c2,d2, a3,b3,c3,d3, a4,b4,c4,d4, a5,b5,c5]
parameters = [x1a, y1a, x1b, y1b, x1c, y1c, x1d, y1d, x1e, y1e, x2a, y2a, x2b, y2b, x2c, y2c, x2d, y2d, x2e, y2e, x3a, y3a, x3b, y3b, x3c, y3c, x3d, y3d, x3e, y3e, x4a, y4a, x4b, y4b, x4c, y4c, x4d, y4d, x4e, y4e, x5a, y5a, x5b, y5b, x5c, y5c, x5d, y5d, x5e, y5e]

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
L3 = [c3;0;d3]
L4 = [c4;0;d4]
L5 = [c5;0;d5]
V = [vx;vy;vz]
W = [wx;wy;wz]

Ld1 = [a1;1;b1]
Ld1x = [0 -b1 1; b1 0 -a1; -1 a1 0]
Ld2 = [a2;1;b2]
Ld2x = [0 -b2 1; b2 0 -a2; -1 a2 0]
Ld3 = [a3;1;b3]
Ld3x = [0 -b3 1; b3 0 -a3; -1 a3 0]
Ld4 = [a4;1;b4]
Ld4x = [0 -b4 1; b4 0 -a4; -1 a4 0]
Ld5 = [a5;1;b5]
Ld5x = [0 -b5 1; b5 0 -a5; -1 a5 0]


u1a = [x1a y1a 1]
u1b = [x1b y1b 1]
u1c = [x1c y1c 1]
u1d = [x1d y1d 1]
u1e = [x1e y1e 1]

u2a = [x2a y2a 1]
u2b = [x2b y2b 1]
u2c = [x2c y2c 1]
u2d = [x2d y2d 1]
u2e = [x2e y2e 1]

u3a = [x3a y3a 1]
u3b = [x3b y3b 1]
u3c = [x3c y3c 1]
u3d = [x3d y3d 1]
u3e = [x3e y3e 1]

u4a = [x4a y4a 1]
u4b = [x4b y4b 1]
u4c = [x4c y4c 1]
u4d = [x4d y4d 1]
u4e = [x4e y4e 1]

u5a = [x5a y5a 1]
u5b = [x5b y5b 1]
u5c = [x5c y5c 1]
u5d = [x5d y5d 1]
u5e = [x5e y5e 1]

eqn1a = u1a * Ld1x * (L1 - y1a*V - y1a^2*W)
eqn1b = u1b * Ld1x * (L1 - y1b*V - y1b^2*W)
eqn1c = u1c * Ld1x * (L1 - y1c*V - y1c^2*W)
eqn1d = u1d * Ld1x * (L1 - y1d*V - y1d^2*W)
eqn1e = u1e * Ld1x * (L1 - y1e*V - y1e^2*W)

eqn2a = u2a * Ld2x * (L2 - y2a*V - y2a^2*W)
eqn2b = u2b * Ld2x * (L2 - y2b*V - y2b^2*W)
eqn2c = u2c * Ld2x * (L2 - y2c*V - y2c^2*W)
eqn2d = u2d * Ld2x * (L2 - y2d*V - y2d^2*W)
eqn2e = u2e * Ld2x * (L2 - y2e*V - y2e^2*W)

eqn3a = u3a * Ld3x * (L3 - y3a*V - y3a^2*W)
eqn3b = u3b * Ld3x * (L3 - y3b*V - y3b^2*W)
eqn3c = u3c * Ld3x * (L3 - y3c*V - y3c^2*W)
eqn3d = u3d * Ld3x * (L3 - y3d*V - y3d^2*W)
eqn3e = u3e * Ld3x * (L3 - y3e*V - y3e^2*W)

eqn4a = u4a * Ld4x * (L4 - y4a*V - y4a^2*W)
eqn4b = u4b * Ld4x * (L4 - y4b*V - y4b^2*W)
eqn4c = u4c * Ld4x * (L4 - y4c*V - y4c^2*W)
eqn4d = u4d * Ld4x * (L4 - y4d*V - y4d^2*W)
eqn4e = u4e * Ld4x * (L4 - y4e*V - y4e^2*W)

eqn5a = u5a * Ld5x * (L5 - y5a*V - y5a^2*W)
eqn5b = u5b * Ld5x * (L5 - y5b*V - y5b^2*W)
eqn5c = u5c * Ld5x * (L5 - y5c*V - y5c^2*W)
eqn5d = u5d * Ld5x * (L5 - y5d*V - y5d^2*W)
eqn5e = u5e * Ld5x * (L5 - y5e*V - y5e^2*W)


eqns = [eqn1a;eqn1b;eqn1c;eqn1d;eqn1e; eqn2a;eqn2b;eqn2c;eqn2d;eqn2e; eqn3a;eqn3b;eqn3c;eqn3d;eqn3e; eqn4a;eqn4b;eqn4c;eqn4d;eqn4e; eqn5a;eqn5b;eqn5c;eqn5d;eqn5e ]
F = System(eqns; variables = variables, parameters = parameters)

#p = rand(length(parameters))
p = randn(ComplexF64,length(parameters))
println(p)
sol = HomotopyContinuation.solve(F;target_parameters=p, show_progress = true)
#sol = monodromy_solve(F;show_progress = true)
print(sol)
