import numpy as np

tau3 = np.array([[1,0],[0,-1]])

def calculate_alpha(t0, mu, Delta):
    return np.array([[2*t0 - mu,Delta],[np.conj(Delta),-2*t0]])

def calculate_beta(t0):
    return -t0 * tau3
