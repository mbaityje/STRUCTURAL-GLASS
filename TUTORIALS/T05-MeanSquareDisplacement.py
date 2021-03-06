#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a trajectory from a .gsd file, and
# calculates the mean square displacement.
#
# To launch a simulation:
# python T5-MeanSquareDisplacement.py trajectory.gsd
#
# For example:
# python T5-MeanSquareDisplacement.py test-output/trajectory.gsd
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import sys
import numpy as np
import gsd.pygsd
import gsd.hoomd
import matplotlib.pyplot as plt

################################################################
#
# FUNCTIONS THAT WOULD GO IN A SEPARATE MODULE
# 
################################################################

def PeriodicSquareDistance(vec_a, vec_b, box_size):
# This function measures the distance between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #substraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return (delta ** 2).sum(axis=-1)

def PeriodicDistance(vec_a, vec_b, box_size):
    return np.sqrt(PeriodicSquareDistance(vec_a, vec_b, box_size))

def CalculateMeanSquareDisplacements(old_pos, new_pos, box_size):
    all_displacements=PeriodicSquareDistance(old_pos, new_pos, box_size)
    return all_displacements.sum()/len(all_displacements)



################################################################
#
# READ ARGUMENTS
# 
################################################################
if len(sys.argv)==2:
    filename = sys.argv[1] #Name of the file with the trajectory
else:
    print("Launch as:")
    print(sys.argv[0]," configuration.gsd")
print(filename)
every_forMemory = 1

################################################################
#
# READ TRAJECTORY
# 
################################################################
with open(filename, 'rb') as flow:
    HoomdFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
    s0=hoomdTraj.read_frame(0) #This is a snapshot of the initial configuration (frame zero)
    Natoms=s0.particles.N
    boxParams=s0.configuration.box
    L=boxParams[0]
    if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
        print('box dimensions are : ', boxParams[0])
        print('and are not supported (only isotropic systems supported)')
        raise SystemExit
    
    Nframes = len(hoomdTraj)
    trajDuration = Nframes/every_forMemory
    print('there are Nframes=', Nframes, 'in the file, but we only use trajDuration=',trajDuration, ' of them.')

    #Now we create a trajectory array, containing the positions of the particles at each t
    trajectory = np.zeros((trajDuration,Natoms,3))
    for time in range(0, trajDuration, 1):
        ## we only look 1 frame every "every_forMemory" frame, to save memory: 
        trajectory[time] = (hoomdTraj[time*every_forMemory].particles.position) # [selectedAtoms]
    HoomdFlow.close()
    
print('Shape of the trajectory array (times, particles, dimensions):', np.shape(trajectory))


################################################################
# 
# CALCULATE MEAN SQUARE DISPLACEMENT
#
################################################################

initialPositions=trajectory[0]
box_size=np.array([L,L,L])
times=np.zeros(trajDuration)
msd=np.zeros(trajDuration)
for t in range(0, trajDuration):
    times[t]=t*every_forMemory
    msd[t]=CalculateMeanSquareDisplacements(initialPositions, trajectory[t], box_size)
#Notice that the previous cycle could have been performed at reading time (but we wanna make the code simple)


################################################################
# 
# FIGURES
#
################################################################

#I merge together the two lists, to produce a nice text output
output_msd=np.column_stack((times, msd))
np.savetxt('./test-output/msd.txt',output_msd,fmt='%g %.14g')
    
#Now I make a plot of the output
plt.loglog(times, msd)
plt.xlabel('t')
plt.ylabel('mean square displacement')
plt.title('Mean square displacement of '+filename)
plt.grid(True)
plt.savefig(filename+"_MSD.png")
plt.show()









