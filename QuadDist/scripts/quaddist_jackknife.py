import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

#Small script for computing average energy and sign


versions = ['A']


dtau = 1/16

curr_dir = os.getcwd()

beta = curr_dir.split("/")[-4][4:]
beta_flt = float(beta.split("p")[0])+float(beta.split("p")[1])/10**len(beta.split("p")[1])
db = int(curr_dir.split("/")[-3][2:])
ele = curr_dir.split("/")[-6]


qmax = 1200

#print(f5_obj.keys())
OPs = ['Z(Q1)%Z', 'Z(Q2)%Z', 'Z(Q3)%Z']
#OPs = ['Q**2']

for OP in OPs:

    therm_arr = []
    acf_arr = []


    filename = '../{}.db{}.beta{}.s10.h5.A'.format(ele,db,beta)
    n_therm = 0

    f5_obj = h5py.File(filename,'r')
  
    dat = f5_obj[OP]
    sign = f5_obj["sign"]

        
    n_samp = dat.shape[0]
    n_proc = dat.shape[1]

    n_meas = n_samp*n_proc
        
        

    qs = dat.shape[2]
    qvals = np.linspace(-qmax+qmax/qs,qmax-qmax/qs,qs)
   
    prob = dat[:,:,:,0].reshape(-1,qs)
    signs = sign[:,:,0].ravel()
    sign_prob = signs[:,np.newaxis]*prob
    final_prob = np.sum(sign_prob, axis=0)/np.sum(signs)

    prob_jk = np.zeros((n_meas,qs))
    for b in range(n_meas):
        prob_jk[b] = (np.sum(sign_prob,axis=0)-sign_prob[b])/(np.sum(signs)-signs[b])

    final_prob_jk = np.sum(prob_jk,axis=0)/n_meas
    err = np.sqrt((n_meas-1)/n_meas)*np.sqrt(np.sum((prob_jk-final_prob_jk)**2,axis=0))


    if OP == 'Z(Q1)%Z':
        op = "Q20"
    if OP == 'Z(Q2)%Z':
        op = "Q2p"
    if OP == 'Z(Q3)%Z':
        op = "Q2m"
        
    np.savetxt(f'{ele}.{op}dist.db{db}.beta{beta}.dat',np.stack([qvals,final_prob,err],axis=1),header="q \t P(q) \t +/-")

        
