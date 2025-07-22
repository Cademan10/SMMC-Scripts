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
        if 'j-proj?' in line:
            v = line.strip('\n').split()[0].split(",")
            print(v)
            max_j = int(float(v[3])/2)
            

    db_val = convert_str_to_ftr(db_val)
    return nt_val, db_val, max_j

def get_hdf5_file():
    files = os.listdir('../.')
    filename = ''

    for f in files:
        if '.h5.A' in f:
            filename = '../' + f
            break
    return filename

########################################################################################################################


HDF5_file = get_hdf5_file()
print('HDF5_file = %s'%HDF5_file)


######################### Set the type of projection you want to do and the number of jackknife bins ######################
projType = 'J'
nbins = 20
#########################################################################################################################

if nbins == 1:
    os.system("mkdir -p {}_proj".format(projType))
    dirname = '{}_proj'.format(projType)
    os.system("rm -r {}/*".format(dirname))
elif nbins > 1:
    os.system("mkdir -p {}_proj_JK".format(projType))
    dirname = '{}_proj_JK'.format(projType)
    os.system("rm -r {}/*".format(dirname))

   
        

for OP in ['M\(tau,{}\)'.format(projType.upper())]:
    
    cmd = 'pfsmmca --decorr 1  --jkbins {} calc_qtau_sproj_corrmtx '.format(nbins) + OP + ' ' + projType.upper() + ' '+ HDF5_file 
    os.system(cmd)

    nt_str, db_ftr, maxj = get_params()


    print('db_ftr = %s'%db_ftr)
    print('nt_str = %s'%nt_str)


    label=''
    if (OP == 'Q\(tau,{}\)'.format(projType.upper())):
        label='qtau'
    elif (OP == 'M\(tau,{}\)'.format(projType.upper())):
        label='mtau'

    if projType == "M":
        loopStart = -maxj
    elif projType == 'J':
        loopStart = 0

    if nbins > 1:
        for j in range(loopStart,maxj+1):
            os.system("mkdir -p {}/{}{}".format(dirname, projType.upper(),j))

    #cmd = 'mv t.qtau.db' + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx t.qtau.db' + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
    #os.system(cmd)

    
 
    for b in range(1,nbins+1):
            
        for j in range(loopStart,maxj+1):

            if nbins==1:
                cmd = 'mv t.{}_{}{}.db'.format(label,projType.lower(),j) + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx {}/t.{}_{}{}.db'.format(dirname,label,projType.lower(),j) + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
                os.system(cmd)
            elif nbins > 1:
                cmd = 'mv t.{}_{}{}_jk{}.db'.format(label,projType.lower(),j,b) + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx {}/{}{}/t.{}_{}{}_jk{}.db'.format(dirname,projType.upper(),j,label,projType.lower(),j,b) + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
                os.system(cmd)




for OP in ['M\(tau\)M']:

    label=''
    if (OP == 'Q\(tau\)Q'):
        label='qtau'
    elif (OP == 'M\(tau\)M'):
        label='mtau'
    
    cmd = 'pfsmmca --decorr 1 calc_qtau_corrmtx ' + OP + ' '+ HDF5_file 
    os.system(cmd)

    
    cmd = 'mv t.mtau.db' + db_ftr + '.nt' + nt_str + '.smmc.response.corrmtx {}/t.mtau_true.db'.format(dirname) + db_ftr + '.nt' + nt_str + '.no_blocking.smmc.response.corrmtx'
    os.system(cmd)


