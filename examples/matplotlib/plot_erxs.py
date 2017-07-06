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

# Calculate T bin width and midpoint
Tbins_midpoint = [0]*(len(Tbins)-1)
Tbins_dT = [0]*(len(Tbins)-1)
for i in range(len(Tbins) - 1):
    Tbins_dT[i] = Tbins[i] - Tbins[i + 1]
    Tbins_midpoint[i] = np.mean([Tbins[i],Tbins[i + 1]])

L = 8
T = len(Tbins) - 1;
G = len(Eg) - 1;
print(L,T,G)

# Read .csv output file
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

# Calculate sum in Legendre polynomials for all g -> t
# function: compute_xs_gt(mu_L, g, t, erxs_coeff_np)
# function: compute_xs_matrix(mu_L, G, T, erxs_coeff_np)
f = plt.figure(1)
m = np.linspace(0,1,101)
mu_L_list = np.ndarray.tolist(m)
#xs_sum = [0]*len(mu_L_list)
g = 0
x = [0]*len(mu_L_list)
t_list = [3,7,14,21]
t_id = 0
for t in t_list:
    for i in range(len(mu_L_list)):
        mu_L = mu_L_list[i]
        xs_matrix = compute_xs_matrix(mu_L, G, T, erxs_coeff_np)
        #print(xs_matrix)
        #print(mu_L,xs_matrix[g][t])
        x[i] = xs_matrix[g][t]

    plt.subplot(4, 1, t_id + 1)
    plt.plot(mu_L_list,x,'bo-')
    plt.xlim(0,1)
    plt.ylim(-0.005,0.025)
    plt.grid()
    t_id += 1

#plt.xlabel('mu_L')
#plt.title('Elastic Recoil XS vs. Recoil Angle Cosine (mu_L) in Lab Frame')
plt.show()
#sys.exit("")
