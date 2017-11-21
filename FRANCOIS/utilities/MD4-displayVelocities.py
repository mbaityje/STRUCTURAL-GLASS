## this computes:
#
import gsd.pygsd
import gsd.hoomd
import numpy as np
import sys
######
## plots: 
import matplotlib.pyplot as plt
######
plt.ion()   ## optional


filename = sys.argv[1]


############################
### read the file (gsd) ####
with open(filename, 'rb') as flow:
    HoomdFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
    #boxParams=(hoomdTraj.read_frame(0)).configuration.box
    #for time in range(0, trajDuration, 1):
    time= 0
    Natoms= len(hoomdTraj[0].particles.position)
    Nframe = min(10000,len(hoomdTraj))
    vels        = np.zeros((Nframe*Natoms,3))
    normsVels   = np.zeros((Nframe*Natoms))
    for i in range(Nframe):
        pos = hoomdTraj[i].particles.position
        vel = hoomdTraj[i].particles.velocity

        normsVel = (vel[:,0]**2+vel[:,1]**2+vel[:,2]**2)**0.5

        vels[i*Natoms:(1+i)*Natoms] = vel[:]
        normsVels[i*Natoms:(1+i)*Natoms] = normsVel[:]


a,b,c=plt.hist(normsVels,50,normed=1)

#a,b = plt.hist(vels,100)


kT=0.43
m=1.0
m_over_2pikT = m/(2*np.pi*kT)
def P_vnorm(v):
    return  m_over_2pikT**1.5 * 4*np.pi*v**2*np.exp(-m*v**2/(2*kT))
plt.figure()
plt.plot(b[:-1],a)    
plt.plot(b, P_vnorm(b), linewidth=5)



