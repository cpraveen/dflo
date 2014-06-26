L   = 10.0  # length of domain
n   = 100   # no. of cells
gam = 1.4

dx = L/n

#core conditions

#Total energy
Et = 3.2e6

#Energy per volume
Ec = Et/(dx*dx)
pc = Ec*(gam-1.0)

#outer conditions
E0 = 1.0e-10
p0 = E0*(gam-1.0)

print "Pressure      in core = ", pc
print "Pressure      outside = ", p0
print "Energy/volume in core = ", Ec
print "Energy/volume outside = ", E0
