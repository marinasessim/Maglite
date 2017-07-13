import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
import sys
from compute import *

# Parameters used
Eg = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]
Tbins = [2097152, 1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0]

#Eg = [2, 1]
#Tbins = [2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

L = 7
T = len(Tbins) - 1;
G = len(Eg) - 1;
print(L,T,G)

# Read .csv output file with cross section coefficients
erxs_coeff_np = np.zeros((L + 1, G, T))
row_id = 0
data = np.loadtxt('case_A_erxs_out.csv', delimiter = ',')
with open('case_A_erxs_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        i = 0
        for g in range(G):
            for l in range(L + 1):
                erxs_coeff_np[l, g, row_id] = data[row_id, i]
                i += 1
        row_id += 1

# Read .csv output file with mu_L_max and mu_L_min
mu_L_min_max = np.zeros((2, G, T))
row_id = 0
data = np.loadtxt('case_A_mu_L_out.csv', delimiter = ',')
with open('case_A_mu_L_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        i = 0
        for g in range(G):
            for j in range(2):
                mu_L_min_max[j, g, row_id] = data[row_id, i]
                i += 1
        row_id += 1

m = np.linspace(0,1,1001)
mu_L_list = np.ndarray.tolist(m)

# We choose beforehand the (g, t) combination we want to plot
g = 3
#t = 8

t_list = [8,10,15,20]
t_id = 0
y = [[], [], [], []]
mu_L = [[], [], [], []]

# Post processor to check if mu_L is in the possible range
t_id = 0
for t in t_list:
    m = np.linspace(mu_L_min_max[1][g][t],mu_L_min_max[0][g][t],1001)
    mu_L_list = np.ndarray.tolist(m)
    delta_mu_L = mu_L_min_max[0][g][t] - mu_L_min_max[1][g][t]
    for i in range(len(mu_L_list)):
        mu_L_scaled = 2.0 * (mu_L_list[i] - mu_L_min_max[1][g][t]) / delta_mu_L - 1
        xs_matrix = compute_xs_matrix(mu_L_scaled, delta_mu_L, L, G, T, erxs_coeff_np)
        y[t_id].append(xs_matrix[g][t])
        mu_L[t_id].append(mu_L_list[i])
    t_id += 1
    print(t_id)

# Two subplots, sharing x axis
f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
ax1.set_title('Sharing X axis')
ax1.plot(mu_L[0], y[0], 'b-')
ax1.grid()
ax2.plot(mu_L[1], y[1], 'k-')
ax2.grid()
ax3.plot(mu_L[2], y[2], 'r-')
ax3.grid()
ax4.plot(mu_L[3], y[3], 'y-')
ax4.grid()
plt.xlim(0,1)

#plt.grid()
plt.show()

sys.exit("")
