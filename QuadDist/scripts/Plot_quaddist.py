import numpy as np
import matplotlib.pyplot as plt
import os


OPs = ["Q20","Q2p","Q2m"]
#OPs = ["Q20"]

curr_dir = os.getcwd()

beta = curr_dir.split("/")[-4][4:]
beta_flt = float(beta.split("p")[0])+float(beta.split("p")[1])/10**len(beta.split("p")[1])
db = int(curr_dir.split("/")[-3][2:])
ele = curr_dir.split("/")[-6]


for op in OPs:
    dat = np.loadtxt(f"{ele}.{op}dist.db{db}.beta{beta}.dat")

    qs = dat[:,0]
    prob = dat[:,1]
    err = dat[:,2]


    plt.errorbar(qs[::5],prob[::5],yerr=err[::5],color="b")
    #   plt.fill_between(qs,prob-err,prob+err,color='b',alpha=0.2)
    plt.ylabel("P(q)")
    plt.xlabel("q")
    plt.title(f"$P(Q_{{2{op.split('2')[1]}}})$ vs $q$ for $^{{{ele[2:]}}}{ele[0:2]}$ at $\\beta = ${beta_flt}")
    plt.savefig(f"{ele}_{op}dist_beta{beta}.png")

    plt.close()
    

