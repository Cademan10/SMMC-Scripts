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
    #path = os.path.dirname(os.getcwd())
    path = os.getcwd()
    vec = path.strip('\n').split('/')
    #idx = vec.index('smmc')
    #nucleus = vec[idx+1]
    #beta_str = vec[idx+2]
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
params= {'text.latex.preamble' : r'\usepackage{amsmath}'}
plt.rcParams.update(params)


nucleus, beta_str, beta_ftr = get_parameters()
betas=[str_to_flt(beta_str)]

datafile = np.ndarray(shape=(2), dtype="S10000")
datafile[0] = 'MtM_decorrelation.txt'
datafile[1] = 'QtQ_decorrelation.txt'



t_vals = []
for line in open(datafile[0]):
    values = line.strip("\n").split()
    it = (int(values[0]))
    t_vals.append(it)

tmin = min(t_vals)
thalf = int(ceil(0.5*max(t_vals)))
tquarter = int(ceil(0.25*max(t_vals)))
diffs = [abs(x - tquarter) for x in t_vals]
tquarter = tquarter + min(diffs)
print("tau min = %s"%tmin)
print("tau half = %s"%thalf)
print("tau max = %s"%(max(t_vals)))
print("tau quarter = %s"%(tquarter))

fig, ax = plt.subplots(1,2, figsize=(10,6))
plt.subplots_adjust(hspace=0, wspace=0)

for j,oper in enumerate(['M','Q']):
    x1 = []
    x2 = []
    y1 = []
    y2 = []
    print("file = %s"%(datafile[j]))
    for line in open (datafile[j]):
            values = line.strip ("\n").split ()
            it = (int(values[0]))
            isamp = (int (values[1]))
            qtau_q = (float (values[2]))
            if(abs(qtau_q) > 1.0):
                continue
            if it == tmin:
                x1.append(int (values[1]))
                y1.append(float (values[2]))
            if it == tquarter:
                x2.append(int (values[1]))
                y2.append(float (values[2]))
    
    b2 = 0.5*betas[0]
    #b4 = 0.25*betas[0]
    b4 = round(tquarter*betas[0]/max(t_vals),3)
    ax[j].plot(x1,y1,".",color='red',label=r'$\tau=0 { \; \rm{MeV} }^{ -1 }$')
    ax[j].plot(x2,y2,".",color='blue',label=r'$\tau=%s { \; \rm{MeV} }^{ -1 }$'%b4)
    ax[j].plot(x2,[0.05 for x in range(len(x2))],"--",color='black')
    ax[j].plot(x2,[-0.05 for x in range(len(x2))],"--",color='black')
ax[1].legend(loc=1,fontsize=10)
    


#ax[0].set_ylim(get_lim(My1,My2))
#ax[1].set_ylim(get_lim(Qy1,Qy2))


#ax[1].ticklabel_format(axis='y', style='sci',scilimits=(0,0),useOffset='False')

ax[1].set_yticklabels([])


ax[0].set_xlabel("Separation")
ax[1].set_xlabel("Separation")

ylbl = r'Autocorrelation'
ax[0].set_ylabel(ylbl)

name, A = decompose_nuclei_A(nucleus)

ax[0].set_title(r'${ {  }^{ %s } }$%s(M1)'%(A,name))
ax[1].set_title(r'${ {  }^{ %s } }$%s(Q2)'%(A,name))


output=nucleus+'_acf_'+beta_str+'_db32.pdf'


#plt.title('Interesting Graph\nCheck it out')
print("File generated %s"%output)
plt.savefig(output)
plt.close()
#plt.show()
