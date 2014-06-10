mach = 0.63
aoa  = 2
rho = 1.0
v   = 1.0
gamma = 1.4
p = 1/(gamma * mach^2)

vx = v * cos(aoa*pi/180);
vy = v * sin(aoa*pi/180);

w_0 = rho
w_1 = vx
w_2 = vy
w_3  = p
