import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
import sys
from compute import *

# Parameters used
#Eg = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]
#Tbins = [2097152, 1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0]

Eg = [2, 1]
Tbins = [2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

L = 7
T = len(Tbins) - 1;
G = len(Eg) - 1;
print(L,T,G)

# Read .csv output file with cross section coefficients
erxs_coeff_np = np.zeros((L + 1, G, T))
row_id = 0
data = np.loadtxt('case_test_shiftedLP_out.csv', delimiter = ',')
with open('case_test_shiftedLP_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        i = 0
        for g in range(G):
            for l in range(L + 1):
                erxs_coeff_np[l, g, row_id] = data[row_id, i]
                i += 1
        row_id += 1

# Read .csv output file with min and max mu_o
# mu_o_max, mu_o_min
mu_L_min_max = np.zeros((2, G, T))
row_id = 0
data = np.loadtxt('mu_L_out.csv', delimiter = ',')
with open('mu_L_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        i = 0
        for g in range(G):
            for j in range(2):
                mu_L_min_max[j, g, row_id] = data[row_id, i]
                i += 1
        row_id += 1

print(mu_L_min_max)

# Calculate sum in Legendre polynomials for all g -> t
# function: compute_xs_gt(mu_L, g, t, erxs_coeff_np)
# function: compute_xs_matrix(mu_L, G, T, erxs_coeff_np)
m = np.linspace(0,1,1001)
mu_L_list = np.ndarray.tolist(m)

# We choose before hand the (g, t) combination we want to plot
g = 0
t_list = [5,10,15,20]
t_id = 1
for t in t_list:
    x = [0]*len(mu_L_list)
    mu_L = [0]*len(mu_L_list)
    print(t)
    for i in range(len(mu_L_list)):
        delta_mu_L = mu_L_min_max[0][g][t] - mu_L_min_max[1][g][t]
        value = 2.0 * (mu_L_list[i] - mu_L_min_max[1][g][t]) / delta_mu_L - 1
        xs_matrix = compute_xs_matrix(value, delta_mu_L, G, T, erxs_coeff_np)
        x[i] = xs_matrix[g][t]

        # Post processor to check if mu_L is in the possible range
        if mu_L_list[i] < mu_L_min_max[0][g][t] and mu_L_list[i] > mu_L_min_max[1][g][t]:
            print(mu_L_list[i],x[i])
            plt.subplot(4, 1, t_id)
            plt.plot(mu_L_list[i],x[i],'bo')
            plt.xlim(0,1)
            plt.grid()
    t_id += 1

plt.show()

#sys.exit("")
