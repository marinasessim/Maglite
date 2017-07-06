import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb
from legendreP import legendreP

# Parameters used
Eg = [2, 1]
Tbins = [2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

mu_L = [0, 0, 0.999931, 0.975302, 0.949285, 0.922535, 0.894986, 0.866561, 0.837172, 0.806713, 0.775058, 0.742053, 0.707511, 0.671193, 0.632794, 0.59191, 0.547983, 0.500214, 0.447373, 0.387391, 0.316228, 0.223447]

G = len(Eg)
T = len(Tbins)
print(G,T)

Tbins_midpoint = [Tbins[i] - 0.05 for i in range(len(Tbins) - 1)]
#print(Tbins_midpoint)

# Read .csv output file and save to erxs_coeff
erxs_coeff = []
with open('case_anlt_out.csv', 'rb') as myfile:
    reader = csv.reader(myfile, delimiter=',')
    for row in reader:
        erxs_coeff.append(row[0:])
print(erxs_coeff)

# Save average cross section in Tbin
erxs_ave = np.resize([],(G - 1, T - 1))
for j in range(G - 1):
    for i in range(T - 1):
        x = float(erxs_coeff[i][j])
        #print x
        y = Tbins[i] - Tbins[i+1]
        #print y
        erxs_ave[j][i] = x/y
        #print erxs_ave

# Plot average cross section in Tbin
f1 = plt.figure(1)
for i in range(len(erxs_ave)):
    plt.plot(Tbins_midpoint[:], erxs_ave[i], 'dr', hold=True) #basex = np.exp(1)
plt.title('Simple case for analytical comparison \n Elastic Recoil XS vs. Recoil Energy')
plt.ylabel('ERXS [barn]')
plt.xlabel('T bins [eV]')
plt.xlim(0,3)
#plt.ylim(0,500)
#plt.legend(['group 1'])
plt.grid()

f2 = plt.figure(2)
for i in range(len(erxs_ave)):
    plt.plot(mu_L[:], erxs_ave[i], 'db', hold=True) #basex = np.exp(1)
plt.title('Simple case for analytical comparison \n Elastic Recoil XS vs. mu_L')
plt.ylabel('ERXS [barn]')
plt.xlabel('mu_L')
plt.xlim(-1,1)
#plt.ylim(0,500)
#plt.legend(['group 1'])
plt.grid()

plt.show()
