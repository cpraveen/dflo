gamma = 1.4

rho_l = 1.0
v_l   = 0.0
p_l   = 1.0
rho_r = 0.125
v_r   = 0.0
p_r   = 0.1


fprintf(1,'Left state\n')
w_0 = rho_l * v_l
w_1 = 0
w_2 = rho_l
w_3 = p_l/(gamma-1.0) + 0.5*rho_l*v_l^2

fprintf(1,'Right state\n')
w_0 = rho_r * v_r
w_1 = 0
w_2 = rho_r
w_3 = p_r/(gamma-1.0) + 0.5*rho_r*v_r^2
