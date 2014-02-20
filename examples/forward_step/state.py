gamma = 1.4

ux = 3
uy = 0
p  = 1
rho= gamma

w_0 = rho * ux
w_1 = rho * uy
w_2 = rho
w_3 = p/(gamma-1) + 0.5*rho*(ux*ux + uy*uy)

print 'set w_0 value = {0:14.5f}'.format(w_0)
print 'set w_1 value = {0:14.5f}'.format(w_1)
print 'set w_2 value = {0:14.5f}'.format(w_2)
print 'set w_3 value = {0:14.5f}'.format(w_3)
