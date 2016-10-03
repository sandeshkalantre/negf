#calculate the DC josephson effect in the tunnelling approximation

import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.constants as const

import helper as hlp
import surface

def calculate_josephson(E, kT, mu, M, g_1, g_2):
    return (const.e)/(const.h) * (np.abs(M)**2) * 4 * np.real(g_1[1,0] * g_2[0,1] - g_2[1,0] * g_1[0,1]) * 1.0/(1.0 + np.exp((E - mu)/kT))

def calculate_I_J_E(t0, mu, kT, M, Delta, phi, N = 200, L = 10):
    # N : number of points in the energy linspace
    # L : number of Delta to calculate away from the gap
    alpha1 = hlp.calculate_alpha(t0, mu, Delta) 
    beta1 = hlp.calculate_beta(t0)

    alpha2 = hlp.calculate_alpha(t0, mu, Delta * np.exp(1j*phi)) 
    beta2 = hlp.calculate_beta(t0)

    E = np.linspace(mu - L*Delta, mu + 2*Delta, N) 
    I_J = np.zeros(N)
    for i in range(N):
        g1 = surface.calculate_surface(E[i], alpha1, beta1) 
        g2 = surface.calculate_surface(E[i], alpha2, beta2) 

        I_J[i] = calculate_josephson(E[i], kT, mu, M, g1, g2) 

    return E, I_J

#energy in eV
t0 = 0.1
Delta = 0.08
kT = 0.002
mu = 0.0
M = 1e-2

phi1 = np.pi/4
phi2 = np.pi/2
phi3 = 5*np.pi/4
phi4 = 3*np.pi/2
E1,I_J1 = calculate_I_J_E(t0,mu,kT,M,Delta,phi1)
E2,I_J2 = calculate_I_J_E(t0,mu,kT,M,Delta,phi2)
E3,I_J3 = calculate_I_J_E(t0,mu,kT,M,Delta,phi3)
E4,I_J4 = calculate_I_J_E(t0,mu,kT,M,Delta,phi4)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(E1/Delta, I_J1,label=r"$\phi = \frac{\pi}{4}$",linewidth=1.5)
ax.plot(E2/Delta, I_J2,label=r"$\phi = \frac{\pi}{2}$",linewidth=1.5)
ax.plot(E3/Delta, I_J3,label=r"$\phi = \frac{5\pi}{4}$",linewidth=1.5)
ax.plot(E4/Delta, I_J4,label=r"$\phi = \frac{3\pi}{2}$",linewidth=1.5)
ax.set_xlabel(r"$\frac{E}{|\Delta|}$", fontsize=24)
ax.set_ylabel(r"$I_J$", fontsize=24)
ax.set_title(r"Energy Spectrum near the superconducting gap")
ax.legend(fontsize=24)
plt.show()


