import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

#Small script for computing average energy and sign
os.chdir("..")

filename = r'MultipleBetas/Ne20.db64.beta2.s10.h5.A'
n_therm = 0

f5_obj = h5py.File(filename,'r')

OP1 = 'Q**2'
OP2 = 'sign'

H = f5_obj[OP1]
sign = f5_obj[OP2]

n_samp = H.shape[0]
n_proc = H.shape[1]

H_acf = []
H_therm = []
num = []

j_max = 100

# Thermalization
for i in range(0,n_samp):
    num.append(i+1)
    H_therm.append(H[i][0][0])

# Average
ha = 0.0
for i in range(n_therm,n_samp):
    ha = ha + H[i][0][0]/(n_samp-n_therm)

# ACF (Average Correlation Function)

for j in range(0,j_max):
    acf = 0.0

    # Calculates the avg. CF of all samples (exlcuding burn-in) seperated by j steps
    for i in range(n_therm,n_samp):
        acf = acf + (H[i][0][0] - ha)*(H[i-j][0][0] - ha)
    H_acf.append(acf/(n_samp-n_therm))
#print(H_acf)
acf0 = H_acf[0]

H_acf[j] = np.array(H_acf[j])/acf0



f = open("energy_acf.dat",'w')
for i in range(0,j_max):
    f.write(str(i) + ' ' + str(H_acf[i]) + '\n')

plt.plot(list(range(0,j_max)),H_acf)
plt.savefig("Figs/avg_corr_func.pdf")
plt.close()

plt.plot(num,H_therm)
plt.savefig("Figs/energy_therm.pdf")
