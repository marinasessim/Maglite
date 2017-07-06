import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
import sys
from compute import *

# Parameters used
Eg = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]
Tbins = [2097152, 1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0]

# Calculate T bin width and midpoint
Tbins_midpoint = [0]*(len(Tbins)-1)
Tbins_dT = [0]*(len(Tbins)-1)
for i in range(len(Tbins) - 1):
    Tbins_dT[i] = Tbins[i] - Tbins[i + 1]
    Tbins_midpoint[i] = np.mean([Tbins[i],Tbins[i + 1]])

L = 9
T = len(Tbins) - 1;
G = len(Eg) - 1;
print(L,T,G)

# Read .csv output file
erxs_coeff_np = np.zeros((L + 1, G, T))
row_id = 0
data = np.loadtxt('test_legendre_10.csv', delimiter = ',')
with open('case_A_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        i = 0
        for g in range(G):
            for l in range(L + 1):
                erxs_coeff_np[l, g, row_id] = data[row_id, i]
                i += 1
        row_id += 1

# Calculate sum in Lengendre polynomials for all g -> t
# function: compute_xs_gt(mu_o, g, t, erxs_coeff_np)
# function: compute_xs_matrix(mu_o, G, T, erxs_coeff_np)
f1 = plt.figure(1)
mu_o_list = [-1 + 0.1 * i for i in range(21)]
#xs_sum = [0]*len(mu_o_list)
g = 1
t = 8
print ('hello world')
for i in range(len(mu_o_list)):
    mu_o = mu_o_list[i]
    xs_matrix = compute_xs_matrix(mu_o, G, T, erxs_coeff_np)

    #print(xs_matrix)
    print(mu_o,xs_matrix[g][t])

    plt.plot(mu_o,xs_matrix[g][t],'ob')
    plt.hold(True)

plt.show()

sys.exit("")

# Save average cross section in Tbin
erxs_ave = np.resize([],(G, T))
for j in range(G):
    for i in range(T):
        x = float(erxs_coeff[i][j])
        y = Tbins[i] - Tbins[i+1]
        erxs_ave[j][i] = x/y

# Plot average cross section in Tbin
f1 = plt.figure(1)
for i in range(len(erxs_ave)):
    plt.semilogx(Tbins_midpoint[:], erxs_ave[i], '-', hold=True) #basex = np.exp(1)
plt.title('Case A \n Hydrogen-1: Elastic Recoil XS Averaged vs. Recoil Energy')
plt.ylabel('ERXS [barn]')
plt.xlabel('T bins [eV]')
#plt.xlim(0,3)
#plt.ylim(0,500)
plt.legend(['group 1', 'group 2', 'group 3', 'group 4', 'group 5', 'group 6', 'group 7'])
plt.grid()

plt.show()
