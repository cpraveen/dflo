import numpy as np


def compute(d):
    for i in range(1,4):
        rhov = np.log(d[i-1,2]/d[i,2])/np.log(2.0)
        rho  = np.log(d[i-1,3]/d[i,3])/np.log(2.0)
        E    = np.log(d[i-1,4]/d[i,4])/np.log(2.0)
        print rhov, rho, E

print "Q1"
d = np.loadtxt("polystate_isoscheme_q1.txt")
compute(d)

print "Q2"
d = np.loadtxt("polystate_isoscheme_q2.txt")
compute(d)
