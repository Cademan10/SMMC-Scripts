import numpy as np
import matplotlib.pyplot as plt
import os


OPs = ["Q20","Q2p"]
#OPs = ["Q20"]

curr_dir = os.getcwd()

ele = curr_dir.split("/")[-4]
beta = curr_dir.split("/")[-2][4:]
beta_flt = float(beta.split("p")[0])+float(beta.split("p")[1])/10**len(beta.split("p")[1])

do_marg_dist = True

for op in OPs:
    dat = np.loadtxt(f"t.{op}.db.0000000.smmc.dat")

    qs = dat[:,0]
    prob = dat[:,1]
    err = dat[:,2]


    plt.errorbar(qs,prob,yerr=err,color="b",label="SMMC")

    if do_marg_dist:
        f = open(f"../../Data/{op.lower()}dist.txt")
        dat = []
        lines = f.readlines()
        for i in range(len(lines)):
            if "beta" in lines[i]:

                if float(lines[i].split()[3]) == beta_flt:
                    data_start = i+2
                    break
        for i in range(data_start,len(lines)):
            if lines[i] == "\n":
                break
            dat.append([float(line) for line in lines[i].split()])
        dat = np.array(dat)

        plt.plot(dat[:,0],dat[:,1],label="Marginal distribution",zorder=1000,color='r',linestyle="--")
        
    #   plt.fill_between(qs,prob-err,prob+err,color='b',alpha=0.2)
    plt.ylabel("P(q)")
    plt.xlabel("q")
    plt.legend()
    plt.title(f"$P(Q_{{2{op.split('2')[1]}}})$ vs $q$ for $^{{{ele[2:]}}}{ele[0:2]}$ at $\\beta = ${beta_flt}")
    plt.savefig(f"{ele}_{op}dist_beta{beta}.png")
    
    plt.close()
    

