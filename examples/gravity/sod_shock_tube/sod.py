import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True

d100q1 = np.loadtxt("sod_100_q1/Density.curve")
d200q1 = np.loadtxt("sod_200_q1/Density.curve")
d100q2 = np.loadtxt("sod_100_q2/Density.curve")
d200q2 = np.loadtxt("sod_200_q2/Density.curve")
fvm    = np.loadtxt("fvm/sod2000.dat")

lw = 2
fs = 18

plt.figure()
plt.plot(d100q1[:,0],d100q1[:,1],'o',linewidth=lw)
plt.plot(d100q2[:,0],d100q2[:,1],'*',linewidth=lw)
plt.plot(fvm[:,0],fvm[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Density",fontsize=fs)
plt.legend(("Q1, 100 cells","Q2, 100 cells","FVM, 2000 cells"))
plt.savefig("sod_rho_100.eps")

plt.figure()
plt.plot(d200q1[:,0],d200q1[:,1],'o',linewidth=lw)
plt.plot(d200q2[:,0],d200q2[:,1],'*',linewidth=lw)
plt.plot(fvm[:,0],fvm[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Density",fontsize=fs)
plt.legend(("Q1, 200 cells","Q2, 200 cells","FVM, 2000 cells"))
plt.savefig("sod_rho_200.eps")
