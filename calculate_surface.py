import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import inv
import numpy.linalg as LA

#test parameters, do not affect the calculation of the surface Greem's function

t0 = 1.0
mu = 0.5
Delta = 2.0 * np.exp(1j * np.pi/2.0)

tau3 = np.array([[1,0],[0,-1]])
alpha = np.array([[2*t0 - mu,Delta],[np.conj(Delta),-2*t0 + mu]])
beta = -t0 * tau3


def calculate_surface(E,alpha,beta,eta = 0.00001,N = 10000,eps = 0.0000000001,brute=False):
    # N : number of iterations
    g0 = alpha + 1j * eta * np.eye(2)
    g = g0
    g_last = g0

    for i in range(N):
        g = inv((E+1j*eta)*np.eye(2) - alpha - np.dot(np.conj(beta).T,np.dot(g,beta)))
        g = 0.5 * (g + g0) 
        if brute == False and LA.norm(g_last - g)/LA.norm(g_last) < eps:
            #print "--------------------------------------------------------------------------------"
            #print "Calculated in",i,"iterations for E = ",E 
            #print "alpha =\n",alpha
            #print "g =\n",g
            #print "\n"
            return g
        else:
            g_last = g
    #print "--------------------------------------------------------------------------------"
    #print "Calculated in",N,"iterations for E = ",E 
    #print "alpha =\n",alpha
    #print "g =\n",g
    #print "\n"
    return g

#print calculate_surface(1.1,alpha,beta)
    

