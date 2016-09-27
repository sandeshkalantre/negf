#calculate the DC josephson effect in the tunnelling approximation

import numpy as np
import matplotlib.pyplot as plt
import time

from helper import *
import calculate_surface

def calculate_josephson(g_1,g_2):
    return 4*np.real(g_1[1,0] * g_2[0,1] - g_2[1,0] * g_1[0,1])
t0 = 10.0
mu = 100.0
Delta = 1.5


def calculate_I_J_phi(t0,mu,Delta):
    st = time.time()

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
    E = np.linspace(mu_1-10*Delta,mu_1+10*Delta,M) #number of points in phase linspace
    N = 100
    phi = np.linspace(0,2*np.pi,N)
    I_J = np.zeros(N)

    for i in range(N):
        print "Calculating for phase difference = ",phi[i]
        for j in range(M):
            Delta_2 = Delta * np.exp(1.0j * phi[i])
            alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_1 + mu_2]])

            g_1 = calculate_surface.calculate_surface(E[j], alpha_1,beta_1)
            g_2 = calculate_surface.calculate_surface(E[j], alpha_2,beta_2)

            I_J[i] += calculate_josephson(g_1, g_2) * 1.0/(1+np.exp((E[j]-mu)/0.001)) 
        #print "\n"

    print "Done for Delta = " + str(Delta) + " in " + str(time.time()-st) + " seconds."
    return phi,I_J

def calculate_I_J_E(t0,mu,Delta,phi=np.pi/2):
    st = time.time()

    t0_1 = t0
    mu_1 = mu
    Delta_1 = Delta * np.exp(1j * 0.0)

    alpha_1 = np.array([[2*t0_1 - mu_1,Delta_1],[np.conj(Delta_1),-2*t0_1 + mu_1]])
    beta_1 = -t0_1 * tau3


    t0_2 = t0
    mu_2 = mu
    Delta_2 = Delta * np.exp(1.0j * phi)

    alpha_2 = np.array([[2*t0_2 - mu_2,Delta_2],[np.conj(Delta_2),-2*t0_2 + mu_2]])
    beta_2 = -t0_2 * tau3
    
    #number of points in Energy linspace
    M = 100
    #since the iteration to calculate surface Green's function blows up near the superconducting gap, the energy values are kept lower than that
    E = np.linspace(mu_1-20*Delta, mu_1+20*Delta, M) #number of points in phase linspace
    I_J = np.zeros(M)

    for i in range(M):
        print "Calculating for Energy = ", E[i]

        g_1 = calculate_surface.calculate_surface(E[i], alpha_1,beta_1,brute=True)
        g_2 = calculate_surface.calculate_surface(E[i], alpha_2,beta_2,brute=True)

        I_J[i] += calculate_josephson(g_1, g_2) * 1.0/(1+np.exp((E[i]-mu)/1.0)) 

    print "Done for Delta = " + str(Delta) + " in " + str(time.time()-st) + " seconds."
    return E,I_J

#phi,I_J = calculate_I_J_phi(t0,mu,0.5)
#phi_1,I_J_1 = calculate_I_J_phi(t0,mu,1.0)
#phi_2,I_J_2 = calculate_I_J_phi(t0,mu,1.5)

E,I_J_3 = calculate_I_J_E(t0,mu,1.5,np.pi/2)
E,I_J_4 = calculate_I_J_E(t0,mu,1.5,np.pi/4)

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(phi*180/np.pi,I_J,label =  r"$|\Delta| = 0.5$", linewidth=1.5)
#ax.plot(phi_1*180/np.pi,I_J_1,label =  r"$|\Delta| = 1.0$", linewidth=1.5)
#ax.plot(phi_2*180/np.pi,I_J_2,label =  r"$|\Delta| = 1.5$", linewidth=1.5)
#ax.set_xlabel(r"$\Delta \phi$ (degrees)", fontsize=16)
#ax.set_ylabel(r"$I_J$ (Arbitrary Units)", fontsize=16)
#ax.set_title("DC Josephson Effect for different values of " + r"$|\Delta|$") 
#ax.text(0.9,0.9,r"$t = $" + str(t0) + "\n" + r"$\mu = $" + str(mu) + "\n" + r"$|\Delta| = $" + str(np.abs(Delta)),ha="center",va="center",transform=ax.transAxes,fontsize=16) 
#plt.legend(prop={'size':24})
#plt.savefig("jos_delta.png",dpi=300)

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.plot((E-mu)/Delta,I_J_3,linewidth = 1.5,label = r"$\phi = \frac{\pi}{2}$" )
ax.plot((E-mu)/Delta,I_J_4,linewidth = 1.5, label = r"$\phi = \frac{\pi}{4}$" )
ax.set_xlabel(r"$\frac{E-\mu}{|\Delta|}$",fontsize=24)
ax.set_ylabel(r"$I_J$ (Arbitrary Units)",fontsize=16)
ax.set_title("Spectrum for different values of " + r"$\phi$") 
plt.legend(prop={'size':24})
fig1.savefig("jos_E.png",dpi=300)

plt.show()


