import numpy as np
import os


# copying pn.sps r0.red and r2.red
os.system('cp ../../../Inter/p.sps .')
os.system('cp ../../../Inter/n.sps .')
os.system('cp ../../../Inter/r0.red .')
os.system('cp ../../../Inter/quad.red .')
os.system('cp ../../../Inter/*.int .')

r0_file = 'r0.red'
proton_file = 'p.sps'
neutron_file = 'n.sps'

nb_proton_states = len(open(proton_file).readlines(  ))
nb_neutron_states = len(open(neutron_file).readlines(  ))
nb_states = nb_proton_states + nb_neutron_states

r0_neutron={}
r0_proton={}
for i,line in enumerate(open(r0_file)):
    values = line.strip ("\n").split ()
    '2nd block in r0 is for neutron, 1st block is for proton.'
    if i >= nb_proton_states:
        for j,val in enumerate(values):
            r0_neutron[(i-nb_proton_states,j)] = val
    for j,val in enumerate(values):
        r0_proton[(i,j)] = val

r0=[]
for i in range(0,nb_states):
    for j in range(0,nb_states):
        if i < nb_proton_states and j < nb_proton_states:
            r0_val = r0_proton[(i,j)]
        elif i >= nb_proton_states and j >= nb_proton_states:
            r0_val = r0_neutron[(i-nb_proton_states,j-nb_proton_states)]
        else:
            r0_val = 0.0
        r0.append(r0_val)

with open('r0.red','w') as out:
    for val in r0:
        out.write(str(val) + '\n')



os.system('cp p.sps second/')
os.system('cp n.sps second/')
os.system('cp r0.red second/')
os.system('cp quad.red second/')
os.system('cp *.int second/')



