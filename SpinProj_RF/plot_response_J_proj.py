import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

#Small script for computing average energy and sign

db = 1/32
beta = 0.625
nt = int(beta/db)


plot_spin_projs = False

#print(f5_obj.keys())
OPs = ['M(tau,J)']
#OPs = ['Q**2']

for OP in OPs:
    
    files = np.array(["J_proj/"+f for f in os.listdir("J_proj/") if str(nt) in f])

    labels = []
    for f in files:
        
        l = f.split(".")[1].split("_")[1]

        if l[0] == 'j':
            labels.append(int(l.split("j")[1]))
        else:
            labels.append(l)
    
    j_srt = np.argsort(np.array([ (0, x) if isinstance(x, (int, float)) else (1, x) for x in labels],dtype=[('x',int),('y',object)]),order=['x','y'])
   
    files = files[j_srt]
    labels = np.array(labels)[j_srt]

    for i in range(len(files)):
        

        dat = open(files[i]).readlines()
        
        qtaus = []
        response = []
        err = []
        
        for lines in dat:
            line = lines.split()


            if line[1] == "Covariance":
                break
            
            if line[0] == "#":
                continue
      
            
            qtaus.append(float(line[0]))
            response.append(float(line[1]))
            err.append(float(line[2]))
        
        if labels[i] == 'true':
            plt.errorbar(qtaus,response,yerr=err,label="J = {}".format(labels[i]),color='k')

        elif labels[i]=='full':
            plt.plot(qtaus,response,label="J = {}".format(labels[i]),zorder=1000,ls='--')
            plt.fill_between(qtaus,np.array(response)+np.array(err),np.array(response)-np.array(err),alpha=0.4)
            
        elif plot_spin_projs:
            plt.plot(qtaus,response,label="J = {}".format(labels[i]),ls="--")
            plt.fill_between(qtaus,np.array(response)+np.array(err),np.array(response)-np.array(err),alpha=0.4)

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4)


    if OP == "M(tau,J)":
        op = "M1"

    if OP == "Q(tau,J)":
        op = "Q2"
    
    plt.title("$R_{{{}}}(\\tau,J)$ for $^{{20}}Ne$ \n $\\beta = {}$".format(op,str(beta)))
    plt.savefig("Figs/{}.png".format(OP),dpi=700,bbox_inches='tight')
    plt.show()
    plt.show()
 #   print(files)
    

        
        
