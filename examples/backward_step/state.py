gamma = 1.4

# Normal shock calaculations can be done here
# http://www.dept.aoe.vt.edu/~devenpor/aoe3114/calc.html 

#Left state
u_l = 4.08 # Obtained from normal shock relation
v_l = 0
p_l  = 30.059 # Obtained from normal shock relation
rho_l = 7.0406 # Obtained from normal shock relation

#Right state
u_r = 0.0
v_r = 0.0
p_r  = 1.0
rho_r = 1.4

#Left state
w_0_l = rho_l * u_l
w_1_l = rho_l * v_l
w_2_l = rho_l
w_3_l = p_l/(gamma-1) + 0.5*rho_l*(u_l*u_l + v_l*v_l)

#Right state
w_0_r = rho_r * u_r
w_1_r = rho_r * v_r
w_2_r = rho_r
w_3_r = p_r/(gamma-1) + 0.5*rho_r*(u_r*u_r + v_r*v_r)

#Print left state
print 'Left state'
print 'w_0 value = {0:14.5f}'.format(w_0_l)
print 'w_1 value = {0:14.5f}'.format(w_1_l)
print 'w_2 value = {0:14.5f}'.format(w_2_l)
print 'w_3 value = {0:14.5f}'.format(w_3_l)

#Print right state
print '\nRight state'
print 'w_0 value = {0:14.5f}'.format(w_0_r)
print 'w_1 value = {0:14.5f}'.format(w_1_r)
print 'w_2 value = {0:14.5f}'.format(w_2_r)
print 'w_3 value = {0:14.5f}'.format(w_3_r)
