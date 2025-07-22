import numpy as np
import math
import os
#from scipy.optimize import curve_fit
import numpy.linalg as nl

#os.system("mkdir -p jk")
def regularize_cov_matrix(data_response,cov_matrix):
    for i in range(len(data_response)):
        cov_matrix[i,i] = pow(data_response[i][2],2)
        #print('[%s] : cov = %s ; data = %s'%(i,cov_matrix[i,i],pow(data_response[i][2],2)))


def extrapolate_matrix_V2(matrix_db32, matrix_db64, C_matrix):
    x1 = 1./64
    x2 = 1./32

    dim = np.shape(matrix_db32)[0]
    matrix_db0 = np.zeros((dim,dim))

    matrix_db32 = matrix_db32.astype(float)
    matrix_db64 = matrix_db64.astype(float)

    for i in range(dim):
        for j in range(dim):
            matrix_db0[i,j] = C_matrix[i,0] * C_matrix[j,0] * matrix_db64[i,j] + C_matrix[i,1] * C_matrix[j,1] * matrix_db32[i,j]
 
 

    return matrix_db0

def extrapolate_response_general(data_db32, data_db64):
    data_db0 = []
    x1 = 1./64
    x2 = 1./32

    dim = len(data_db32)
    C_matrix = np.zeros((dim,2))

    dX = x2-x1

    for i in range(len(data_db32)):
        E = float(data_db32[i][0])
        y1 = float(data_db64[i][1])
        y2 = float(data_db32[i][1]) 
        dY = y2 - y1

        errdb64 = float(data_db64[i][2])
        errdb32 = float(data_db32[i][2]) 
        derr = errdb32 - errdb64

        sigma1 = errdb64
        sigma2 = errdb32

        Sxx = (pow(x1,2))/(pow(sigma1,2))
        Sxx += (pow(x2,2))/(pow(sigma2,2))
        S = (1./pow(sigma1,2)) + (1./pow(sigma2,2))
        Sx = (x1/pow(sigma1,2)) + (x2/pow(sigma2,2))

        Delta = S*Sxx - Sx*Sx
        iDelta = 1./Delta

        C1 = iDelta*(Sxx - Sx*x1)/pow(sigma1,2)
        C2 = iDelta*(Sxx - Sx*x2)/pow(sigma2,2)
        C_matrix[i,0] = C1
        C_matrix[i,1] = C2

        y0 = C1*y1 + C2*y2 

        # Derivate Formula
        # a + b*x
        # b = (y2 - y1)/dx
        # a = y1 - (y2 - y1)*x1/dx
        # da = (1 + x1/dx)dy1 - (x1/dx)dy2
        #err0 = (1 + x1/dX)*errdb64 - (x1/dX)*errdb32

        err0 = Sxx/ Delta

        data_db0.append( (E, y0, err0) )

    return data_db0, C_matrix

def extrapolate_response(data_db32, data_db64):
    data_db0 = []
    x1 = 1./64
    x2 = 1./32

    dim = len(data_db32)
    C_matrix = np.zeros((dim,2))

    dX = x2-x1

    for i in range(len(data_db32)):
        E = float(data_db32[i][0])
        y1 = float(data_db64[i][1])
        y2 = float(data_db32[i][1]) 
        dY = y2 - y1

        errdb64 = float(data_db64[i][2])
        errdb32 = float(data_db32[i][2]) 
        derr = errdb32 - errdb64

        sigma1 = errdb64
        sigma2 = errdb32

        by = dY/dX

        C1 = (x2*x2 - x1*x2)/pow(x1-x2,2)
        C2 = (x1*x1 - x1*x2)/pow(x1-x2,2)
        C_matrix[i,0] = C1
        C_matrix[i,1] = C2

        y0 = C1*y1 + C2*y2 

        # Derivate Formula
        # a + b*x
        # b = (y2 - y1)/dx
        # a = y1 - (y2 - y1)*x1/dx
        # da = (1 + x1/dx)dy1 - (x1/dx)dy2
        #err0 = (1 + x1/dX)*errdb64 - (x1/dX)*errdb32

        Sxx = (pow(x1,2))/(pow(sigma1,2))
        Sxx += (pow(x2,2))/(pow(sigma2,2))
        
        S = (1./pow(sigma1,2)) + (1./pow(sigma2,2))

        Sx = (x1/pow(sigma1,2)) + (x2/pow(sigma2,2))

        Delta = S*Sxx - Sx*Sx
        err0 = Sxx/ Delta

        data_db0.append( (E, y0, np.sqrt(err0)) )

    return data_db0, C_matrix


def get_data(filename):
    lines = [line.strip('\n').rstrip() for line in open(filename)]
    idx_cov = lines.index('# Covariance matrix')
    idx_tau = lines.index('# tau   Q(tau)   +/-')
    
    header = lines[0:idx_tau+1]
    response = lines[idx_tau+1:idx_cov]
    cov = lines[(idx_cov+1):]

    data_response = [x.split() for x in response]
    data_cov = [x.split() for x in cov]

    return data_response, data_cov, header
    

####################################################################################################

do_jk = True
projType = "J".upper()


dirName = f"{projType}_proj"

if do_jk:
    dirName += "_JK"

os.system("mkdir -p %s"%dirName)
os.system("rm -r {}/*".format(dirName))

dirs = []

if not do_jk:
   
    for d in os.listdir(f"../dB32/second/sf/{dirName}"):

        dirs.append(" ")
    
    

if do_jk:
    for d in os.listdir(f"../dB32/second/sf/{dirName}"):
        #        os.system("mkdir -p {}/jk_{}".format(dirName, b))
        if "{}".format(projType) in d:
            os.system("mkdir {}/{}".format(dirName,d))
            dirs.append(d)

dbetas = ['0156250','0312500']
op = 'mtau'

####################################################################################################

# Extrapolate covariances for jackknife samples
for d in dirs:
    dname = dirName

    if do_jk:
        dname = dirName + "/"+d

    list_files = os.listdir(f'../dB32/second/sf/{dname}/')
    file_dbeta32 = []


    for f in list_files:
        if (op in f) and len(f.split(" ")) == 1 and len(f.split("#")) == 1:
            file_dbeta32.append(f'../dB32/second/sf/{dname}/' + f)
    list_files = os.listdir(f'../dB64/second/sf/{dname}/')
    file_dbeta64 = []
   # print(file_dbeta32)
    for f in list_files:
        if (op in f) and len(f.split(" ")) == 1 and len(f.split("#")) == 1:
            file_dbeta64.append(f'../dB64/second/sf/{dname}/' + f)


    if do_jk:
        jks_32 = [int(f.split("/")[-1].split("_")[2].split(".")[0][2:]) for f in file_dbeta32]

        file_dbeta32 = np.array(file_dbeta32)[np.argsort(jks_32)]

        jks_64 = [int(f.split("/")[-1].split("_")[2].split(".")[0][2:]) for f in file_dbeta64]
        file_dbeta64 = np.array(file_dbeta64)[np.argsort(jks_64)]

        jks = np.array(jks_32)[np.argsort(jks_32)]

    else:
        js_32 = [f.split("/")[-1].split("_")[1].split(".")[0] for f in file_dbeta32]
        js_32 = np.array(js_32)[np.argsort(js_32)]
        js_64 = [f.split("/")[-1].split("_")[1].split(".")[0] for f in file_dbeta64]
        js_64 = np.array(js_32)[np.argsort(js_64)]
        js = np.array(js_32)[np.argsort(js_32)]
       
       # file_dbeta32 = np.array(file_dbeta32)[np.argsort(jks_32)]
       # jks_64 = [int(f.split("/")[-1].split("_")[2].split(".")[0][2:]) for f in file_dbeta64]
       # file_dbeta64 = np.array(file_dbeta64)[np.argsort(jks_64)]
       # jks = np.array(jks_32)[np.argsort(jks_32)]
    
    for i in range(len(file_dbeta32)):
#    for i in range(1):
        data_response_db32, data_cov_db32, header = get_data(file_dbeta32[i])
  
        data_response_db64, data_cov_db64, _ = get_data(file_dbeta64[i])

        # Rearrange 1D array into 2D array matrix
        data_cov_db32 = np.array(data_cov_db32)
        data_cov_db64 = np.array(data_cov_db64)



        data_response_db0, C_matrix = extrapolate_response(data_response_db32, data_response_db64)
        cov_matrix_db0 = extrapolate_matrix_V2(data_cov_db32, data_cov_db64,C_matrix)



        Es0, vec0 = nl.eig(cov_matrix_db0)
        for k,e in enumerate(Es0):
            if e < 0:
                print('Problem! Negative eigenvalue in extrapolated covariance matrix for %s.'%op)
                print('Eigenvalue [%s] = %s'%(k,e))
            if (e > 0) and (e < 1e-6):
                print('Zero eigenvalue [%s] = %s for %s'%(k,e,op))
                #print('Eigenvalue [%s] = %s'%(i,e))
        print(dname)
        if do_jk:
            file_db0 = '{}/t.'.format(dname)+op+'.extrapolated.jk{}'.format(jks[i])+'.smmc.response.corrmtx'
        else:
            file_db0 = '{}/t.'.format(dname)+op+'.extrapolated.{}'.format(js[i])+'.smmc.response.corrmtx'
        header[3] = '# db -> 0 extrapolation'
  
        with open(file_db0,'w') as out:
            for line in header:
                out.write(line + '\n')
            for i in range(len(data_response_db0)):
                E = data_response_db0[i][0]
                res = data_response_db0[i][1]
                err = data_response_db0[i][2]
                out.write(str(E) + ' ' + str(res) + ' ' + str(err) + '\n' )
            out.write('# Covariance matrix'+'\n')
            for i in range(len(data_response_db0)):
                cov_line = ''
                for j in range(len(data_response_db0)):
                    cov_line += str(cov_matrix_db0[i,j]) + ' '

                cov_line += '\n'
                out.write(cov_line)

            data_cov_db32 = data_cov_db32.astype(float)
            data_cov_db64 = data_cov_db64.astype(float)
        out.close()


if do_jk:

    # Extrapolate covariances of full dataset
    list_files = os.listdir('../dB32/second/sf/{}/'.format(dirName))
    file_dbeta32 = ''
    for f in list_files:
        if (op in f) and ("true" in f):
            file_dbeta32 = '../dB32/second/sf/{}/'.format(dirName) + f

    list_files = os.listdir('../dB64/second/sf/{}/'.format(dirName))
    file_dbeta64 = ''
    for f in list_files:
        if (op in f) and ("true" in f):
            file_dbeta64 = '../dB64/second/sf/{}/'.format(dirName) + f

    data_response_db32, data_cov_db32, header = get_data(file_dbeta32)
    data_response_db64, data_cov_db64, _ = get_data(file_dbeta64)

    # Rearrange 1D array into 2D array matrix
    data_cov_db32 = np.array(data_cov_db32)
    data_cov_db64 = np.array(data_cov_db64)

    data_response_db0, C_matrix = extrapolate_response(data_response_db32, data_response_db64)
    cov_matrix_db0 = extrapolate_matrix_V2(data_cov_db32, data_cov_db64,C_matrix)

    Es0, vec0 = nl.eig(cov_matrix_db0)
    for i,e in enumerate(Es0):
        if e < 0:
            print('Problem! Negative eigenvalue in extrapolated covariance matrix for %s.'%op)
            print('Eigenvalue [%s] = %s'%(i,e))
        if (e > 0) and (e < 1e-6):
            print('Zero eigenvalue [%s] = %s for %s'%(i,e,op))
            #print('Eigenvalue [%s] = %s'%(i,e))

    file_db0 = '{}/t.'.format(dirName)+op+'.extrapolated.true.smmc.response.corrmtx'

    header[3] = '# db -> 0 extrapolation'

    with open(file_db0,'w') as out:
        for line in header:
                out.write(line + '\n')
        for i in range(len(data_response_db0)):
            E = data_response_db0[i][0]
            res = data_response_db0[i][1]
            err = data_response_db0[i][2]
            out.write(str(E) + ' ' + str(res) + ' ' + str(err) + '\n' )
        out.write('# Covariance matrix'+'\n')
        for i in range(len(data_response_db0)):
            cov_line = ''
            for j in range(len(data_response_db0)):
                cov_line += str(cov_matrix_db0[i,j]) + ' '

            cov_line += '\n'
            out.write(cov_line)

    data_cov_db32 = data_cov_db32.astype(float)
    data_cov_db64 = data_cov_db64.astype(float)


    os.system('cp {}/t.'.format(dirName)+op+'.extrapolated.true.smmc.response.corrmtx .')
