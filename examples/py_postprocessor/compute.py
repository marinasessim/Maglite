import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
from legendreP import *
import sys
import math

# Computes the summation of the expansion in Legendre for given t and g
def compute_xs_gt(mu_L_scaled, delta_mu_L, L, g, t, erxs_coeff_np):
    #L = erxs_coeff_np.shape[0] - 1
    G = erxs_coeff_np.shape[1]
    T = erxs_coeff_np.shape[2]
    s = 0
    for l in range(L + 1):
        s += (2 * l + 1) / delta_mu_L * erxs_coeff_np[l,g,t] * legendreP(l, mu_L_scaled)
    return s

# Computes the cross section matrix for all T's and G's
def compute_xs_matrix(mu_L_scaled, delta_mu_L, L, G, T, erxs_coeff_np):
    temp = np.zeros((G, T))
    for t in range(T):
        for g in range(G):
            temp[g,t] = compute_xs_gt(mu_L_scaled, delta_mu_L, L, g, t, erxs_coeff_np)
    return temp
