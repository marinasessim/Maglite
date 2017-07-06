import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
from legendreP import legendreP

# Parameters used
Eg = [2, 1]
Tbins = [10.0 - 0.1 * i for i in range(101)]
Tbins_midpoint = [Tbins[i] - 0.05 for i in range(100)]
#print(Tbins)
print(Tbins_midpoint)

T = len(Tbins)
G = len(Eg)
#print(T,G)

# Read .csv output file
erxs_coeff = []
with open('case_B.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        erxs_coeff.append(row[0:])
#print(erxs_coeff)

# Save average cross section in Tbin
erxs_ave = np.resize([],(G - 1, T - 1))
for j in range(G - 1):
    for i in range(T - 1):
        x = float(erxs_coeff[j][i])
        y = np.mean([Tbins[i],Tbins[i+1]])
        erxs_ave[j][i] = x/y

# Plot cross section vs. mid point in T bins
f1 = plt.figure(1)
for i in range(len(erxs_coeff)):
    plt.plot(Tbins_midpoint[:], erxs_coeff[i], 'd', hold=True) #basex = np.exp(1)
plt.title('Case B \n Carbon-12: Elastic Recoil XS vs. Recoil Energy')
plt.ylabel('ERXS [barn]')
plt.xlabel('T bins [eV]')
plt.xlim(0,3)
#plt.ylim(0,500)
#plt.legend(['group 1'])
plt.grid()

# Plot average cross section in Tbin
f2 = plt.figure(2)
for i in range(len(erxs_ave)):
    plt.plot(Tbins_midpoint[:], erxs_ave[i], 'dr', hold=True) #basex = np.exp(1)
plt.title('Case B \n Carbon-12: Elastic Recoil XS Averaged vs. Recoil Energy')
plt.ylabel('ERXS [barn]')
plt.xlabel('T bins [eV]')
plt.xlim(0,3)
#plt.ylim(0,500)
#plt.legend(['group 1'])
plt.grid()

plt.show()
