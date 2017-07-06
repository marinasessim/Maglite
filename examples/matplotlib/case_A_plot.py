import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
from legendreP import legendreP

# Parameters used
Eg = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]
Tbins = [2097152, 1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0]

# Calculate T bin width and midpoint
Tbins_midpoint = [0]*(len(Tbins)-1)
Tbins_dT = [0]*(len(Tbins)-1)

for i in range(len(Tbins) - 1):
    Tbins_dT[i] = Tbins[i] - Tbins[i + 1]
    Tbins_midpoint[i] = np.mean([Tbins[i],Tbins[i + 1]])

T = len(Tbins)
G = len(Eg)

# Read .csv output file
erxs_coeff = []
with open('case_A_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        erxs_coeff.append(row[0:])

# Save average cross section in Tbin
erxs_ave = np.resize([],(G - 1, T - 1))
for j in range(G - 1):
    for i in range(T - 1):
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
