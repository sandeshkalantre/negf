#calculate the DC josephson effect in the tunnelling approximation

import numpy as np
import matplotlib.pyplot as plt

from helper import *
import calculate_surface

def calculate_josephson(g_1,g_2):
    return 4*np.real(g_1[1,0] * g_2[0,1] - g_2[1,0] * g_1[0,1])

t0_1 = 4.0
mu_1 = 10.0
Delta_1 = 2.0 * np.exp(1j * 0.0)

alpha_1 = np.array([[2*t0_1 - mu_1,Delta_1],[np.conj(Delta_1),-2*t0_1 + mu_1]])
beta_1 = -t0_1 * tau3


t0_2 = 4.0
mu_2 = 10.0
Delta_2 = 2.0 * np.exp(1.0j * 0.0)

alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_2 + mu_2]])
beta_2 = -t0_2 * tau3

#number of points in Energy linspace
M = 50
#since the iteration to calculate surface Green's function blows up near the superconducting gap, the energy values are kept lower than that
E = np.linspace(mu_1-10,mu_1+2,M)
#number of points in phase linspace
N = 100
phi = np.linspace(0,2*np.pi,N)
I_J = np.zeros(N)

for i in range(N):
    print "Calculating for phase difference = ",phi[i]
    for j in range(M):
        Delta_2 = 2.0 * np.exp(1.0j * phi[i])
        alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_1 + mu_2]])

        g_1 = calculate_surface.calculate_surface(E[j], alpha_1,beta_1)
        g_2 = calculate_surface.calculate_surface(E[j], alpha_2,beta_2)

        I_J[i] += calculate_josephson(g_1, g_2)
    print "\n"

plt.plot(phi,I_J)
plt.show()

