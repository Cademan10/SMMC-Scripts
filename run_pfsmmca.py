import numpy as np
import os

def convert_str_to_ftr(label_str):

    #label_ftr = ''
    print("String=",label_str)
    if float(label_str) < 1:
        vec = label_str.split('.')
        integer = vec[0]
        decimals = vec[1]
        for i in range(len(decimals),7):
            decimals += '0'
        label_ftr = integer + '.' + decimals
    else:
        decimals = '0000000'
        label_ftr += '.' + decimals

    if float(label_ftr) < 1.:
        label_ftr = label_ftr[1:]

    return label_ftr


def get_params():
    files = os.listdir('../.')
    filename = ''

    for f in files:
        if ('.in' in f) and ('db' in f):
            filename = '../' + f

            if ("~" in f):
                continue

            break

    print('Input file : %s'%filename)
    for line in open(filename):
        if 'dbeta' in line:
            v = line.strip('\n').split()[0].split(",")
            print(v)
            nt_val = v[0]
            db_val = v[1]

            print("dB=",db_val)

    db_val = convert_str_to_ftr(db_val)
    return nt_val, db_val

def get_hdf5_file():
    files = os.listdir('../.')
    filename = ''

    for f in files:
        if '.h5.A' in f:
            filename = '../' + f
            break
    return filename

####################


HDF5_file = get_hdf5_file()
print('HDF5_file = %s'%HDF5_file)

for OP in ['Q\(tau\)Q','M\(tau\)M']:

    label=''
    if (OP == 'Q(tau)Q'):
    	label='qtau'
    elif (OP == 'M(tau)M'):
    	label='mtau'
    
    cmd = 'pfsmmca --decorr 1 calc_qtau_corrmtx ' + OP + ' ' + HDF5_file 
    os.system(cmd)

nt_str, db_ftr = get_params()

print('db_ftr = %s'%db_ftr)
print('nt_str = %s'%nt_str)

cmd = 'mv t.qtau.db' + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx t.qtau.db' + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
os.system(cmd)                                                                           
cmd = 'mv t.mtau.db' + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx t.mtau.db' + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
os.system(cmd)

cmd = 'cp *.smmc.response.corrmtx ../../../sf/'
os.system(cmd)
