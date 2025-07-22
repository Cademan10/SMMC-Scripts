import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import os

os.chdir(r"/Users/cader/Desktop/Yale Research/Research/Ne20_test/SMMC")

font = { 'size'   : 14, 'family' : 'serif'}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : '\\usepackage{amsmath}'}
plt.rcParams.update(params)


SMMC_DATA = np.loadtxt(r"/Users/cader/Desktop/Yale Research/Research/Fe56/smmc/Data/Fe56_data.dat")
HF_SHELL_DATA = np.loadtxt("/Users/cader/Desktop/Yale Research/HF SHELL/examples/Ne20/HF_SHELL_DATA.dat")

Ts_HF = 1/HF_SHELL_DATA[0:,0]
Q2_HF = HF_SHELL_DATA[0:,1]
Q20_HF = HF_SHELL_DATA[0:,2]
H_HF = HF_SHELL_DATA[0:,3]


betas = SMMC_DATA[0:,0]

Q2 = SMMC_DATA[0:,5]
Q2_err = SMMC_DATA[0:,6]

H = SMMC_DATA[0:,1]
H_err = SMMC_DATA[0:,2]

Cv = SMMC_DATA[0:,3]
Cv_err = SMMC_DATA[0:,4]



Ts = 1/np.array(betas[1:])


plt.plot(Ts_HF,Q2_HF,'ks--',label="HF-SHELL")


plt.errorbar(Ts,Q2[1:],yerr=Q2_err[1:],c='b',marker='o',label="SMMC",capsize=2)
plt.legend()
plt.tick_params(axis='both', direction='in')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)

plt.xlabel(r'$T$ (MeV)')
plt.ylabel(r'$\langle Q^2 \rangle_T$ (fm$^4$)')
plt.title('$\\langle Q^2 \\rangle_T$ vs $T$ for $^{20}$Ne')
plt.savefig('Figs/Ne20.Q2.vs.T.png',dpi=800)

plt.close()

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

r0 = 1.2  #fm
A = 20
Chi = 3*r0**2*A**(5/3)/np.sqrt(5*np.pi)

beta = np.sqrt(Q2[1:])/Chi
#beta_HF=np.sqrt(Q2_HF)/Chi
#beta_HF = np.sqrt(5*np.pi)/(3*A**(5/3)*r0**2)*Q20_HF
beta_HF = np.abs(Q20_HF)/Chi
beta_HF_2 = np.sqrt(Q2_HF)/Chi

plt.errorbar(Ts,beta,yerr=Q2_err[1:]/(2*Chi*np.sqrt(Q2[1:])),c='b',marker='o',label="SMMC",capsize=2)
plt.plot(Ts_HF,beta_HF,'ks--',label="HF-SHELL")


plt.legend()
plt.tick_params(axis='both', direction='in')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)

plt.xlabel(r'$T$ (MeV)')
plt.ylabel(r'$\beta$')
plt.title('$\\beta$ vs $T$ for $^{20}$Ne')
plt.savefig('Figs/Ne20.betaDeformation.png',dpi=800)

plt.close()
#plt.show()

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

#plt.errorbar(Ts,H,yerr=H_err,c='b',marker='o',label="SMMC",capsize=2)
gs_inds = np.where(np.array(betas)>=2)
E_gs = np.mean(np.array(H)[gs_inds])
Cv_gs = np.mean(np.array(Cv)[gs_inds])
beta_gs = np.mean(np.array(betas)[gs_inds])

E_gs_err = np.sqrt(np.mean(np.array(H_err)[gs_inds]**2))
Cv_gs_err = np.sqrt(np.mean(np.array(Cv_err)[gs_inds]**2))
print('Ground state energy = ',E_gs,"+/-",E_gs_err)

f = open("Scripts/lvl2.in")
lines = f.readlines()
lines[7] = 'E0 ='+str(round(E_gs,6))+'\t # SMMC g.s. energy'
f.close()


plt.plot(np.array(betas)[gs_inds],E_gs*np.ones(len(np.array(betas)[gs_inds])),"r--",zorder=1000)
plt.errorbar(betas,H,yerr=H_err,c='b',marker='o',label="SMMC",capsize=2)
plt.errorbar(beta_gs,E_gs,yerr=E_gs_err,marker="*",c="r",markersize=10,capsize=2,zorder=2000)
plt.plot(1/Ts_HF,H_HF,'k--',label='HF-SHELL')
plt.tick_params(axis='both', direction='in')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.legend()
plt.xlabel(r'$\beta$ (MeV)$^{-1}$')
plt.ylabel(r'$\langle H \rangle_T$ (MeV)')
plt.title('$\\langle H \\rangle_T$ vs $\\beta$ for $^{20}$Ne')
plt.savefig('Figs/Ne20.Energy.png',dpi=800)

plt.close()

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
plt.plot(np.array(betas)[gs_inds],Cv_gs*np.ones(len(np.array(betas)[gs_inds])),"r--",zorder=1000)
plt.errorbar(betas,Cv,yerr=Cv_err,c='b',marker='o',label="SMMC",capsize=2)
plt.errorbar(beta_gs,Cv_gs,yerr=Cv_gs_err,marker="*",c="r",markersize=10,capsize=2,zorder=2000)
plt.legend()
plt.tick_params(axis='both', direction='in')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)

plt.xlabel(r'$\beta$ (MeV)$^{-1}$')
plt.ylabel(r'$C_v$')
plt.title('$C_v$ vs $\\beta$ for $^{20}$Ne')
plt.savefig('Figs/Ne20.Cv.png',dpi=800)

plt.close()