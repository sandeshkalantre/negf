import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import inv
import numpy.linalg as LA
import pdb

#helper function
def calc_g(E, alpha, beta, g, eta = 1e-8):
    return inv((E + 1j*eta)*np.eye(2) - alpha - np.dot(np.conj(beta).T, np.dot(g, beta)))

def calculate_surface(E, alpha, beta, eta = 1e-8, N = 10000, eps = 1e-6, brute = False):
    # N : number of iterations
    # eta : 0+ added to get the correct retarded Green's function
    # eps : tolerance in the differnce in the norm of the final result divided by norm of the final result
    # brute : if brute is True, N iterations are made, otherwise the iterations continue until eps is reached
   
    #initial guess
    g0 = inv(alpha) 

    g = g0
    g_last = g0

    for i in range(N):
        g = calc_g(E, alpha, beta, g) 
        g = 0.5*(g_last + g)
        if brute == False and LA.norm(g_last - g)/LA.norm(g) < eps:
            return g
        else:
            g_last = g
    return g

    

