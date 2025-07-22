import numpy as np
from math import ceil
import matplotlib.pyplot as plt
import matplotlib
import argparse
import os

def decompose_nuclei_A(nuc):
    'Take, for example Nd150, and return [Nd, 150]'
    idx = 0
    for i,c in enumerate(nuc):
        if c.isdigit():
            idx = i
            break

    return nuc[:idx], nuc[idx:]

def convert_str_to_ftr(label_str):

    label_ftr = ''
    for param in ['beta','eta']:
        if param in label_str:
            label_ftr = label_str.replace(param,'')
            if 'p' in label_ftr:
                vec = label_ftr.split('p')
                integer = vec[0]
                decimals = vec[1]
                for i in range(len(decimals),7):
                    decimals += '0'
                label_ftr = integer + '.' + decimals
            else:
                decimals = '0000000'
                label_ftr += '.' + decimals
            break

    if float(label_ftr) < 1.:
        label_ftr = label_ftr[1:]

    return label_ftr


def get_parameters():
#    path = os.path.dirname(os.getcwd())
    path = os.getcwd()
    vec = path.strip('\n').split('/')
   
#    idx = vec.index('smmc')
    nucleus = vec[-5]
    beta_str = vec[-3]

    beta_ftr = convert_str_to_ftr(beta_str)
    return nucleus, beta_str, beta_ftr

def str_to_flt(b_str):
    b = 0.
    if 'p' in b_str:
        b = float(beta_str.replace('beta','').replace('p','.'))
    else:
        b = float(beta_str.replace('beta',''))
    return b

font = {'size'   : 18, 'family' : 'serif'}
matplotlib.rc('font', **font )
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : '\\usepackage{amsmath}'}
plt.rcParams.update(params)



nucleus, beta_str, beta_ftr = get_parameters()
betas=[str_to_flt(beta_str)]

datafile = np.ndarray(shape=(2), dtype="S10000")
datafile[0] = 'MtM_thermalization.txt'
datafile[1] = 'QtQ_thermalization.txt'


MaxSamp = 800

t_vals = []
for line in open(datafile[0]):
    values = line.strip("\n").split()
    it = (int(values[0]))
    t_vals.append(it)

tmin = min(t_vals)
thalf = int(ceil(0.5*max(t_vals)))
tquarter = int(ceil(0.25*max(t_vals)))
print("tau min = %s"%tmin)
print("tau half = %s"%thalf)
print("tau max = %s"%(max(t_vals)))


mod=5
lim=MaxSamp


fig, ax = plt.subplots(1,2, figsize=(10,6))
plt.subplots_adjust(hspace=0, wspace=0)

for j,oper in enumerate(['M','Q']):
    x1 = []
    x2 = []
    y1 = []
    y2 = []
    erry1=[]
    erry2=[]
    print("file = %s"%(datafile[j]))
    for line in open (datafile[j]):
            values = line.strip ("\n").split ()
            it = (int(values[0]))
            isamp = (int (values[1]))
            if it == tmin and isamp%mod < 1e-6 and isamp <= lim:
                x1.append(int (values[1]))
                y1.append(float (values[2]))
                erry1.append(float(values[3]))
            if it == tquarter and isamp%mod < 1e-6 and isamp <= lim:
                x2.append(int (values[1]))
                y2.append(float (values[2]))
                erry2.append(float(values[3]))
    
    b2 = 0.5*betas[0]
    b4 = 0.25*betas[0]
    ax[j].errorbar(x1,y1,yerr=erry1,linestyle='dotted',color='red',label=r'$\tau=0 { \; \rm{MeV} }^{ -1 }$')
    ax[j].errorbar(x2,y2,yerr=erry2,linestyle='dotted',color='blue',label=r'$\tau=%s { \; \rm{MeV} }^{ -1 }$'%b4)
    



#ax[0].set_ylim(get_lim(My1,My2))
#ax[1].set_ylim(get_lim(Qy1,Qy2))


#for i in range(len(betas)):
#    #ax[i][1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
#    #ax[i][1].set_yscale('log')
#    ax[i][1].yaxis.tick_right()
#    ax[i][1].text(200,200000,r'${\beta=%s \; {\rm{MeV}}^{-1}}$'%(betas[i]),fontsize=15)
#    #ax[i][1].set_yticklabels([])
#    ax[i,j].legend(loc=4,fontsize=10)

ax[1].yaxis.tick_right()

ax[0].legend(loc=4,fontsize=10)

#ax[1].ticklabel_format(axis='y', style='sci',scilimits=(0,0),useOffset='False')


ax[0].set_xlabel("Samples")
ax[1].set_xlabel("Samples")

ylbl = r'$\rm{\mathcal{O}}(\tau)\rm{\mathcal{O}}(0)$'
ax[0].set_ylabel(ylbl)

name, A = decompose_nuclei_A(nucleus)

ax[0].set_title(r'${ {  }^{ %s } }$%s $(\mathcal{O}=M1)$'%(A,name))
ax[1].set_title(r'${ {  }^{ %s } }$%s $(\mathcal{O}=Q2)$'%(A,name))


output=nucleus+'_therm_'+beta_str+'_db32.pdf'


#plt.title('Interesting Graph\nCheck it out')
print("File generated %s"%output)
plt.savefig(output)
