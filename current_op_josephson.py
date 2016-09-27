#calculate the DC josephson effect in the tunnelling approximation

import numpy as np
import matplotlib.pyplot as plt

from helper import *
import calculate_surface

t0 = 10.0
mu = 10.0
Delta = 1.0

t0_1 = t0
mu_1 = mu
Delta_1 = Delta * np.exp(1j * 0.0)

alpha_1 = np.array([[2*t0_1 - mu_1,Delta_1],[np.conj(Delta_1),-2*t0_1 + mu_1]])
beta_1 = -t0_1 * tau3


t0_2 = t0
mu_2 = mu
Delta_2 = Delta * np.exp(1.0j * 0.0)

alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_2 + mu_2]])
beta_2 = -t0_2 * tau3

#number of points in Energy linspace
M = 50
#since the iteration to calculate surface Green's function blows up near the superconducting gap, the energy values are kept lower than that
E = np.linspace(mu_1-10*Delta,mu_1+10*Delta,M)
#number of points in phase linspace
N = 100
phi = np.linspace(0,2*np.pi,N)
I_J = np.zeros(N,dtype=np.complex128)


for i in range(N):
    print "Calculating for phase difference = ",phi[i]
    for j in range(M):
        Delta_2 = Delta * np.exp(1.0j * phi[i])
        alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_1 + mu_2]])

        g_1 = calculate_surface.calculate_surface(E[j], alpha_1,beta_1)
        g_2 = calculate_surface.calculate_surface(E[j], alpha_2,beta_2)

        a1 = 1*np.dot(g_1,np.dot(beta_1,np.dot((g_1 - np.conj(g_1).T),np.dot(beta_1,np.conj(g_1).T))))
        a2 = 1*np.dot(g_2,np.dot(beta_2,np.dot((g_2 - np.conj(g_2).T),np.dot(beta_2,np.conj(g_2).T))))

        I_op = np.dot(tau3,np.dot(g_2,np.dot(tau3,a1))) + np.dot(tau3,np.dot(a2,np.dot(tau3,np.conj(g_1).T))) - np.dot(tau3,np.dot(a1,np.dot(tau3,np.conj(g_2).T))) - np.dot(tau3,np.dot(g_1,np.dot(tau3,a2))) 

        I_J[i] += I_op[0,0] 
    print "\n"

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(phi,I_J)
ax.set_xlabel(r"$\Delta \phi$", fontsize=16)
ax.set_ylabel(r"$I_J$ (Arbitary Units)", fontsize=16)
ax.set_title("Josephson Effect at zero bias")
ax.text(0.9,0.9,r"$t = $" + str(t0) + "\n" + r"$\mu = $" + str(mu) + "\n" + r"$|\Delta| = $" + str(np.abs(Delta)),ha="center",va="center",transform=ax.transAxes,fontsize=16) 
plt.show()

