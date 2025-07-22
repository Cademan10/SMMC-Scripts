import numpy as np
import matplotlib.pyplot as plt

margDist = False
projType = 1

with open("Ne20.db32.beta1p0.s10.t.E","r") as f:
    lines = f.readlines()

    data = []


    for i in range(len(lines)):
        if "projection type {}".format(projType) in lines[i]:
            data_start = i+3
            break
 
    for line in lines[data_start:]:
        
        if "Q-Projected" in line:
            data = data[:-1]
            break
    
        data.append(line.strip().split())
#    data = f.readlines()[167:268]

#print(np.array(data))

data = np.array([[float(value) for value in line] for line in data])

plt.plot(data[:,0],data[:,1])
plt.fill_between(data[:,0],data[:,1]+data[:,2],data[:,1]-data[:,2],alpha=0.5,color='orange')


exact = np.arange(-8*np.sqrt(6),8*np.sqrt(6)+1,np.sqrt(6))


y_maxes = []
for eig in exact:
    y_max = (data[:,1][np.where(np.abs(data[:,0] - eig) == np.min(np.abs(data[:,0] - eig)))]/(np.max(data[:,1])+0.1))[0]
    plt.scatter(eig,y_max*(np.max(data[:,1])+.1),c='r')
    y_maxes.append(y_max)
   # plt.axvline(x=eig,color='r',ymax=y_max,linestyle='--',zorder=0)



if margDist:
    
    margData = np.loadtxt("QuadDist/scripts/q20dist.txt")
    

    
    plt.plot(margData[:,0],margData[:,1])

    print("Oth moment of marg. dist:",np.sum((margData[:-1,1]+margData[1:,1])/2*np.diff(margData[:,0])))
    print("1st moment of marg. dist:",np.sum((margData[:-1,0]+margData[1:,0])/2*(margData[:-1,1]+margData[1:,1])/2*np.diff(margData[:,0])))
    print("2nd moment of marg. dist:",np.sum(((margData[:-1,0]+margData[1:,0])/2)**2*(margData[:-1,1]+margData[1:,1])/2*np.diff(margData[:,0])))
    print("3rd moment of marg. dist:",np.sum(((margData[:-1,0]+margData[1:,0])/2)**3*(margData[:-1,1]+margData[1:,1])/2*np.diff(margData[:,0])))
    print("4th moment of marg. dist:",np.sum(((margData[:-1,0]+margData[1:,0])/2)**4*(margData[:-1,1]+margData[1:,1])/2*np.diff(margData[:,0])))
   # plt.fill_between(data[:,0],data[:,1]+data[:,2],data[:,1]-data[:,2],alpha=0.5,color='orange')




plt.show()
dist_data = np.stack([exact,y_maxes],axis=1)
np.savetxt('Q_dist.dat',dist_data)

print("Oth moment:",np.sum((data[:-1,1]+data[1:,1])/2*np.diff(data[:,0])))
print("1st moment:",np.sum((data[:-1,0]+data[1:,0])/2*(data[:-1,1]+data[1:,1])/2*np.diff(data[:,0])))
print("2nd moment:",np.sum(((data[:-1,0]+data[1:,0])/2)**2*(data[:-1,1]+data[1:,1])/2*np.diff(data[:,0])))
print("3rd moment:",np.sum(((data[:-1,0]+data[1:,0])/2)**3*(data[:-1,1]+data[1:,1])/2*np.diff(data[:,0])))
print("4th moment:",np.sum(((data[:-1,0]+data[1:,0])/2)**4*(data[:-1,1]+data[1:,1])/2*np.diff(data[:,0])))
#print(np.sqrt(np.sum((data[:-1,0]*data[:-1,2]*np.diff(data[:,0]))**2)))



#### Plots the Q-projected energy
with open("Ne20.db32.beta1p0.s10.t.B","r") as f:
    lines = f.readlines()

    data = []


    for i in range(len(lines)):
        if "Q-Projected energy <H>_q for projection type {}".format(projType) in lines[i]:
            data_start = i+3
            break
 
    for line in lines[data_start:]:
        
        if "projection type" in line:
            data = data[:-1]
            break
    
        data.append(line.strip().split())
#    data = f.readlines()[167:268]

#print(np.array(data))

data = np.array([[float(value) for value in line] for line in data])

plt.plot(data[:,0],data[:,1])
plt.fill_between(data[:,0],data[:,1]+data[:,2],data[:,1]-data[:,2],alpha=0.5,color='orange')
plt.ylim(-40,-20)


plt.show()
dist_data = np.stack([exact,y_maxes],axis=1)
np.savetxt('H(Q)_dist.dat',dist_data)


