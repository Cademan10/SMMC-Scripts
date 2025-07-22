import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import os

font = { 'size'   : 14, 'family' : 'serif'}
matplotlib.rc('font', **font)												
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : '\\usepackage{amsmath}'}
plt.rcParams.update(params)

os.system("python lvl2.py lvl2.in")


leveldensity_data = np.loadtxt("../Data/leveldensity.dat")
HF_leveldensity_data = np.loadtxt("/Users/cader/Desktop/Yale Research/HF SHELL/examples/Ne20/HF_LevelDensity.dat")

Ex = leveldensity_data[0:,1]
inds = np.where(Ex>0)

Ex = Ex[inds]

rho = np.exp(leveldensity_data[0:,2])[inds]
print(rho)
rho_err = (np.exp(leveldensity_data[0:,2])*leveldensity_data[0:,3])[inds]




Ex_HF = HF_leveldensity_data[2:,0]
rho_HF = HF_leveldensity_data[2:,2]

plt.errorbar(Ex,rho,yerr=rho_err,label="SMMC",c='b',marker='s',markerfacecolor='none',ls='none')
plt.plot(Ex_HF,rho_HF,'k--',label='HF-SHELL')
plt.xlabel("$E_x$")
plt.ylabel("$\\rho(E)$")
plt.legend()
plt.title('$\\rho$ vs $E_x$ for $^{20}$Ne')
plt.yscale('log') 
plt.savefig('../Figs/Ne20.leveldensity.png',dpi=800)

plt.close()
