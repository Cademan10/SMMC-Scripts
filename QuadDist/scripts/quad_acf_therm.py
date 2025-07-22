import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

#Small script for computing average energy and sign


#versions = ['A','B']
versions = ['E'] 


#print(f5_obj.keys())
OPs = ['Q22-**2','Q22+**2','Q20**2']
#OPs = ['Q**2']

for OP in OPs:

    therm_arr = []
    acf_arr = []
    for l in versions:

        filename = '../Ne20.db32.beta1p0.s10.h5.{}'.format(l)
        n_therm = 0

        f5_obj = h5py.File(filename,'r')
        
        dat = f5_obj[OP]
        print(dat)
        n_samp = dat.shape[0]
        n_proc = dat.shape[1]

        dat_acf = []
        dat_therm = []
        num = []

        j_max = 200

    # Thermalization
        for i in range(0,n_samp):
            num.append(i+1)
            dat_therm.append(dat[i][0][0])


    # Average

        ha = 0.0
        for i in range(n_therm,n_samp):
            ha = ha + dat[i][0][0]/(n_samp-n_therm)

    # ACF (Average Correlation Function)

        for j in range(0,j_max):
            acf = 0.0

        # Calculates the avg. CF of all samples (exlcuding burn-in) seperated by j steps
            for i in range(n_therm,n_samp):
                acf = acf + (dat[i][0][0] - ha)*(dat[i-j][0][0] - ha)
            dat_acf.append(acf/(n_samp-n_therm))
                #print(H_acf)
       
        acf0 = dat_acf[0]

        dat_acf = np.array(dat_acf)/acf0

        acf_arr.append(dat_acf)
        therm_arr.append(dat_therm)

        f = open("{}_v{}_acf.dat".format(OP,l),'w')
        for i in range(0,j_max):
            f.write(str(i) + ' ' + str(dat_acf[i]) + '\n')
    

    for i in range(len(versions)):
        if versions[i] == 'A':
            plt.plot(list(range(0,j_max)),acf_arr[i],label='No angle averaging')

        if versions[i] == 'B':
            plt.plot(list(range(0,j_max)),acf_arr[i],label='Angle averaging')

        else:
            plt.plot(list(range(0,j_max)),acf_arr[i])

    plt.legend()
   
    plt.xlabel('Sample')
    plt.ylabel('Decorrelation')
    plt.title('Decorrelation of {}'.format(OP))
    plt.savefig("{}_acf.png".format(OP))
    plt.close()

    for i in range(len(versions)):
        if versions[i] == 'A':
            plt.plot(num,therm_arr[i],label='No angle averaging')
        if versions[i]== 'B':
            plt.plot(num,therm_arr[i],label='Angle averaging')
        else:
             plt.plot(num,therm_arr[i])
    plt.legend()
    plt.xlabel('Sample')
    plt.ylabel('Observable')
    plt.title('Thermalization of {}'.format(OP))
    plt.savefig("{}_therm.png".format(OP))
    plt.close()
