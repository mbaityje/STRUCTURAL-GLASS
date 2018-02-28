#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# calculates the pair correlation function.
#
# To launch a simulation:
# python T4-PairCorrelation.py configuration.gsd
#
# For example:
# python T4-PairCorrelation.py sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import kobandersen #Here I have put the Kob-Andersen parameters
import matplotlib.pyplot as plt

################################################################
#
# FUNCTIONS THAT WOULD GO IN A SEPARATE MODULE
# 
################################################################

def PeriodicDistance(vec_a, vec_b, box_size):
# This function measures the distance between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #substraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return np.sqrt((delta ** 2).sum(axis=-1))


def CalculateRelativeDistances(positions):
#Takes a vector of n positions.
#Returns a vector of n(n-1)/2 distances
    distances=np.zeros(Natoms*(Natoms-1)/2,dtype=np.double)
    pos=0
    for i in range(Natoms-1):
        distances[pos:pos+Natoms-1-i]+=PeriodicDistance(positions[i],positions[i+1:],np.array([Lx,Ly,Lz]))
        pos+=Natoms-i-1
    return distances

def CalculatePairCorrelationFunction(distances, Natoms, dr=0.1,rMax=2, number_density=1.0):
#Calculate pair correlation function from a clean list of the distances
#between particles (distances between (a,b), (b,a) should appear only once)
#Input:
# -A list of Natoms(Natoms-1)/2 positive distances
#Output:
# -values of g(r): gofr
#
    nbins=int(rMax/dr)+1
    assert(dr>0 and dr<L/2)
    assert(rMax<=L/2)
    print("In the function, Natoms= ",Natoms)
    print(len(distances))
    gofrRaw = np.zeros(nbins)

    for r in distances:
        ibin=int(r/dr)
        if ibin<nbins:
            gofrRaw[int(r/dr)]+=1
    gofrRaw /= (Natoms-1)*0.5 #Normalize by the number of bonds per particle (usually others divide by Natoms, but because they are counting each bond twice)
            
    rvalues = np.arange(0,rMax,dr)
    rvalues[0] = 1   ## correct the fisrt bin, to avoid /0
    gofr = gofrRaw / ( (4*np.pi*rvalues**2*dr) ) #Normalize by the volume of the shell
    gofr /= number_density #So that without structure g(r) goes to 1

    return gofr


################################################################
#
# USUAL STUFF THAT WE SHOWED IN THE PREVIOUS TUTORIALS
# 
################################################################
#Start hoomd
hoomd.context.initialize('--notice-level=0')

#Read arguments in a more simple way
if len(sys.argv)==2:
    filename = sys.argv[1]
else:
    print("Launch as:")
    print(sys.argv[0]," configuration.gsd")
print(filename)

#Read configuration
system = hoomd.init.read_gsd(filename=filename)
Natoms = len(system.particles)
Box = system.box
Lx=Box.Lx
Ly=Box.Ly
Lz=Box.Lz
if(Lx==Ly==Lz):
    L=Lx
assert(Box.get_volume()==Lx*Ly*Lz)
rho=Natoms/(Box.get_volume())

################################################################
# 
# CALCULATE PAIR CORRELATION FUNCTION
#
################################################################

#ACCESS PARTICLE POSITIONS
#To access the array of the particle positions I need to take a
#snapshot (don't ask why it can't be accessed from the system
#directly, but it can't). Anyhow, it's just one more line of code.
snapshot=system.take_snapshot(dtype='double')
positions=np.array(snapshot.particles.position) #Now 'positions' is a NatomsxDIM vector,
                                      #storing all the particles' positions

distances=CalculateRelativeDistances(positions) #A list of the Natoms*(Natoms-1)/2 relative distances

#I set a couple of parameters for the binning of the g(r)
rMax=L/2
dr=0.03

gofr=CalculatePairCorrelationFunction(distances, Natoms, dr=dr, rMax=rMax, number_density=rho)
rvalues = np.arange(0,rMax,dr)


################################################################
# 
# FIGURES
#
################################################################
#I merge together the two lists, to produce a nice text output
output_gofr=np.column_stack((rvalues, gofr))
np.savetxt('./test-output/gofr.txt',output_gofr,fmt='%g %.14g')


#I plot the function I found
plt.plot(rvalues, gofr)
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Pair Correlation function of '+filename)
plt.grid(True)
plt.savefig(filename+"_PAIR.png")
plt.show()











