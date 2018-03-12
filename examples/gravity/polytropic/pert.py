import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.autolayout'] = True


d100q1 = np.loadtxt("pert1_100_q1/Pressure.curve")
d1000q1 = np.loadtxt("pert1_1000_q1/Pressure.curve")

d50q2 = np.loadtxt("pert1_50_q2/Pressure.curve")
d500q2 = np.loadtxt("pert1_500_q2/Pressure.curve")

lw = 2
fs = 18

def initc(x):
    return 1.0e-2 * np.exp(-100*(x-0.5)**2)

plt.figure()
plt.plot(d100q1[:,0],initc(d100q1[:,0]),'--',linewidth=lw)
plt.plot(d100q1[:,0],d100q1[:,1],'o',linewidth=lw)
plt.plot(d1000q1[:,0],d1000q1[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Pressure perturbation",fontsize=fs)
plt.legend(("Initial","100 cells","1000 cells"))
plt.savefig("poly_pert1_q1.eps")

plt.figure()
plt.plot(d100q1[:,0],initc(d100q1[:,0]),'--',linewidth=lw)
#plt.plot(d100q2[:,0],d100q2[:,1]-np.exp(-d100q2[:,0]),'o',linewidth=lw)
plt.plot(d50q2[:,0],d50q2[:,1],'o',linewidth=lw)
plt.plot(d500q2[:,0],d500q2[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Pressure perturbation",fontsize=fs)
plt.legend(("Initial","50 cells","500 cells"))
plt.savefig("poly_pert1_q2.eps")

# pert2
# Pressure is already perturbation pressure

d100q1 = np.loadtxt("pert2_100_q1/Pressure.curve")
d1000q1 = np.loadtxt("pert2_1000_q1/Pressure.curve")

d50q2 = np.loadtxt("pert2_50_q2/Pressure.curve")
d500q2 = np.loadtxt("pert2_500_q2/Pressure.curve")

def initc2(x):
    return 1.0e-4 * np.exp(-100*(x-0.5)**2)


plt.figure()
plt.plot(d100q1[:,0],initc2(d100q1[:,0]),'--',linewidth=lw)
plt.plot(d100q1[:,0],d100q1[:,1],'o',linewidth=lw)
plt.plot(d1000q1[:,0],d1000q1[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Pressure perturbation",fontsize=fs)
plt.legend(("Initial","100 cells","1000 cells"))
plt.savefig("poly_pert2_q1.eps")

plt.figure()
plt.plot(d100q1[:,0],initc2(d100q1[:,0]),'--',linewidth=lw)
#plt.plot(d100q2[:,0],d100q2[:,1],'o',linewidth=lw)
plt.plot(d50q2[:,0],d50q2[:,1],'o',linewidth=lw)
plt.plot(d500q2[:,0],d500q2[:,1],'-',linewidth=lw)
plt.xlabel("x",fontsize=fs)
plt.ylabel("Pressure perturbation",fontsize=fs)
plt.legend(("Initial","50 cells","500 cells"))
plt.savefig("poly_pert2_q2.eps")
