import numpy as np
import os


hdf5_file = '../' + [f for f in os.listdir('../.') if ('h5.A' in f)][0]

cmd = 'smmca_field therm_obs \'M(tau)M\' '+hdf5_file+' ; mv therm_output.txt MtM_thermalization.txt '
os.system(cmd)

cmd = 'smmca_field therm_obs \'Q(tau)Q\' '+hdf5_file+' ; mv therm_output.txt QtQ_thermalization.txt '
os.system(cmd)

