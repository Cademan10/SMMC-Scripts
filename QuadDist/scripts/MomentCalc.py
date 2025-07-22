import numpy as np
import matplotlib.pyplot as plt
import os


OPs = ["Q20","Q2p","Q2m"]
#OPs = ["Q20"]

curr_dir = os.getcwd()

ele = curr_dir.split("/")[-4]
beta = curr_dir.split("/")[-2][4:]
beta_flt = float(beta.split("p")[0])+float(beta.split("p")[1])/10**len(beta.split("p")[1])
n_mom = 4

for op in OPs:
    dat = np.loadtxt(f"t.{op}.db.0000000.smmc.dat")

    qs = dat[:,0]
    prob = dat[:,1]
  
    err = dat[:,2]

    
    moms = [np.sum((prob[0:-1]+prob[1:])/2*((qs[:-1]+qs[1:])/2)**n*np.diff(qs)) for n in range(0,n_mom+1)]
    moms_err = [np.sqrt(np.sum(((err[0:-1]+err[1:])/2*((qs[:-1]+qs[1:])/2)**n*np.diff(qs)))**2) for n in range(0,n_mom+1)]

    np.savetxt(f"{ele}_{op}dist_moments_beta{beta}.dat",np.stack([range(0,n_mom+1),moms,moms_err],axis=1),header="N \t Nth Moment \t +/-")
