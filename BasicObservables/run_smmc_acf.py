import numpy as np
import os


hdf5_file = '../' + [f for f in os.listdir('../.') if ('h5.A' in f)][0]

cmd = 'smmca_field acf_obs --ntherm 200 \'M(tau)M\' '+hdf5_file+' ; mv decorrelation_output.txt MtM_decorrelation.txt '
os.system(cmd)

cmd = 'smmca_field acf_obs --ntherm 200 \'Q(tau)Q\' '+hdf5_file+' ; mv decorrelation_output.txt QtQ_decorrelation.txt '
os.system(cmd)

