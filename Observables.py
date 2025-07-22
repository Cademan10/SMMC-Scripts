import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import os


os.chdir("..")


def getBeta(j):
        os.chdir("dB64")
        
        fName = r'Fe56.db64.beta'+str(j)+'.s10.t.A'
        f = open(fName,"r")

        lines = f.readlines()
        beta = float(lines[2].split()[1])

        if j=="0p0":
                return 0
        else:
                return beta

### obs = 19 for Q2, obs = 8 for H, obs=68 for Cv
def Extrapolation(obs,beta,betaAvg):
        os.chdir("beta"+beta)
        os.chdir("dB64/second")
        fName = r'U238.db64.beta'+str(beta)+'.s10.t.A'
        f = open(fName,"r")

        lines = f.readlines()
#        beta_64 = float(lines[2].split()[1])
       # print(beta_64)
       # print(lines[obs])
        O_err_64 = float(lines[obs].split()[3])
        O_64 = float(lines[obs].split()[1])

	#print(lines[obs].split())

        f.close()

        os.chdir("../../dB32/second")
        
        fName = r'U238.db32.beta'+str(beta)+'.s10.t.A'
        f = open(fName,"r")

        lines = f.readlines()
        beta_32 = float(lines[2].split()[1])
	
        O_err_32 = float(lines[obs].split()[3])
        O_32 = float(lines[obs].split()[1])

	
        f.close()
        os.chdir("../../../")
        
        
        if beta_32<=betaAvg:
       	        return 2*O_64-O_32,np.sqrt(4*O_err_64**2+O_err_32**2)
        else:
	        return (O_64+O_32)/2,np.sqrt(O_err_64**2/4+O_err_32**2/4)

def deformationBeta(beta):
    Q2,Q2_err = Extrapolation(19,beta,50)
        
    r0 = 1.2  #fm
    A = 20
    Chi = 3*r0**2*A**(5/3)/np.sqrt(5*np.pi)

    return (np.sqrt(Q2)/Chi,Q2_err/(2*Chi*np.sqrt(Q2)))

betas = ['40p0','2p0625']
for beta in betas:
    print(deformationBeta(beta))


        
betas = []
Q2 = []
Q2_err = []

H = []
H_err = []

Cv = []
Cv_err = []



dirs = os.listdir(home)
dirs = [dir_i for dir_i in dirs if dir_i[0]=="B"]

for i in range(len(dirs)):
     print("#######")
     dir_i = dirs[i]
     b = dir_i[4:]
     print(b)
     os.chdir(dir_i)
    
     Q2_i,Q2_err_i = Extrapolation(19,b,2)
     H_i,H_err_i = Extrapolation(8,b,2)

     Cv_i,Cv_err_i = Extrapolation(68,b,2)

     beta = getBeta(b)


     betas.append(beta)
     Q2.append(Q2_i)
     Q2_err.append(Q2_err_i)

     H.append(H_i)
     H_err.append(H_err_i)

     Cv.append(Cv_i)
     Cv_err.append(Cv_err_i)

     os.chdir(home)

print(len(betas),len(H),len(Cv),len(Q2))
print("#####")
print(len(H_err),len(Cv_err),len(Q2_err))

sort_inds = np.array(betas).argsort()
betas = np.array(betas)[sort_inds]
H = np.array(H)[sort_inds]
H_err = np.array(H_err)[sort_inds]
Cv = np.array(Cv)[sort_inds]
Cv_err = np.array(Cv_err)[sort_inds]
Q2 = np.array(Q2)[sort_inds]
Q2_err = np.array(Q2_err)[sort_inds]

data = np.stack([betas,H,H_err,Cv,Cv_err,Q2,Q2_err],axis=1)
np.savetxt("Data/Fe56_data.dat",data, delimiter="\t", header="Beta\tH\tH_err\tCv\tCv_err\tQ**2\tQ**2_err")

levelDensity_data = np.stack([betas,H,H_err,Cv,Cv_err],axis=1)

np.savetxt("Data/Fe56_leveldensity_data.dat", levelDensity_data, delimiter="\t", header="Beta\tH\tH_err\tCv\tCv_err")
