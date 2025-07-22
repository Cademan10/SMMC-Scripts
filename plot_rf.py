import numpy as np
import matplotlib.pyplot as plt
import os


files = [f for f in os.listdir('.') if ('db.0000000.no_blocking.smmc.response.corrmtx' in f)]

beta = str(int([f for f in os.listdir('.') if ('t.mtau.db.0156250.' in f)][0].split(".")[4][2:])*0.0156250)




for file in files:
	OP = file.split(".")[1][0]
	#print(OP)
	#beta = file.split(".")[2][4:]+"."+file.split(".")[3][0:3]


	with open(file) as f:
		dat = []
		lines = f.readlines()

		for line in lines:
			#print(line.strip().split())
			if line.strip().split()[0]=="#":
				if line.strip().split()[1]=="Covariance":
					break
				else:
					continue
			else:
				dat.append(line)

	#print(dat)
	dat = np.loadtxt(dat)
	

	tau = dat[0:,0]
	strength = dat[0:,1]
	unc = dat[0:,2]
	
	plt.errorbar(tau,strength,yerr=unc,marker="o",c="b",markersize=3)
	plt.xlabel("$\\tau \\,\\, ($MeV$^{-1})$")

	if OP=="m":
		plt.ylabel("$R_{M1}\\,\\, (\\mu_N^2\\,$MeV$^{-1})$")

	if OP=="q":
		plt.ylabel("$R_{Q2}$")

	plt.title(f"$R_{{{OP}}}$ vs $\\tau$ for $^{{56}}$Fe \n $\\beta$ = "+beta)
	
	plt.savefig("rf_{op}_beta{w}p{p}.png".format(op=OP,w=beta.split(".")[0],p=beta.split(".")[1]),dpi=700)

	plt.close()