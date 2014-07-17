import math

gam = 1.4

# density
rho_l = 8.0
rho_r = 1.4

# velocity
theta = 30.0*math.pi/180.0
u_l = 8.25 * math.cos(theta)
v_l =-8.25 * math.sin(theta)

u_r = 0.0
v_r = 0.0

# pressure
p_l = 116.5
p_r = 1.0

# left state
print "Left state:"
print "set w_0 value = ", rho_l * u_l
print "set w_1 value = ", rho_l * v_l
print "set w_2 value = ", rho_l
print "set w_3 value = ", p_l/(gam-1.0) + 0.5*rho_l*(u_l**2 + v_l**2)

# right state
print "Right state:"
print "set w_0 value = ", rho_r * u_r
print "set w_1 value = ", rho_r * v_l
print "set w_2 value = ", rho_r
print "set w_3 value = ", p_r/(gam-1.0) + 0.5*rho_r*(u_r**2 + v_r**2)
