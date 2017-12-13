#!/usr/bin/python
import numpy as np
from numba import jit
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy import optimize as opt

#-------------------------#
###########################
#        PROFILING        #
###########################
#-------------------------#
import time
#
# Measures time elapsed in function. 
# Example of usage:
# def stampa(a,b,c):
#    print("dada")
#    print(a)
#    print(b)
#    print(c)
# profile(stampa,('ddd','sss','eee'))
#
def profile(func, fargs):
    start = time.time()
    output=func(*fargs)
    end = time.time()
    print("Time elapsed: ",end - start," s")
    return output



#-------------------------#
###########################
#        DISTANCES        #
###########################
#-------------------------#

#------------------------------------
# -- Between particles v_i and v_j --
#------------------------------------
# ParticleSqDist between particles v_i and v_j:
# d_{ij}^2 = \sum_\alpha^D (v_{i,\alpha} - v_{j,\alpha} )^2
# ParticleDist:
# d_{ij} = \sqrt( \sum_\alpha^D (v_{i,\alpha} - v_{j,\alpha} )^2 )
#
######################################################################
# Euclidean distance between two positions (generally two 3D vectors
# indicating particles). Can be feeded 2 lists of positions. In that
# case returns a list of distances.
def ParticleDistPBC(particle_a, particle_b, box_size):
    return np.sqrt(ParticleSqDistPBC(particle_a,particle_b,box_size))

######################################################################
# Euclidean square distance between two positions (generally two 3D
# vectors). Can be feeded 2 lists of positions. In that case returns a
# list of square distances.
def ParticleSqDistPBC(particle_a, particle_b, box_size):
    delta = np.abs(particle_a - particle_b) #usual distance
    delta = np.where(delta > box_size-delta, delta - box_size, delta) #apply PBC
    return np.square(delta).sum(axis=-1)



#-------------------------------------------------------
# -- Distance between configurations at time t and t' --
#-------------------------------------------------------
# ConfDistPBC (slow) and PeriodicDistance (1000 times faster) both calculate:
# Delta = \sum_{t,t'} d_{t,t'}

######################################################################
# Calculates the distance between the same configuration at two different times.
# The operation is not vectorized, so it is slow. I used it for debug
# because I tend to mess up with even the simplest vectorized operations.
# Rather use PeriodicDistance().
def ConfDistPBC(positions_a, positions_b, box_size):
    sum=0
    for i in range(len(positions_a)):
        sum+=ParticleDistPBC(positions_a[i],positions_b[i],box_size)
    return sum

######################################################################
# Calculates the distance between the same configuration at two different times.
# The operation is vectorized, so use this one instead of ConfDistPBC()
def PeriodicDistance(vec_a, vec_b, box_size):
    delta = np.abs(vec_a - vec_b) #subtraction and abs are done component by component
    delta = np.where(delta > box_size-delta, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return np.sqrt(np.square(delta).sum(axis=1)).sum()



#--------------------------------------------------------------
# -- Square Distance between configurations at time t and t' --
#--------------------------------------------------------------
# ConfSqDistPBC (slow) and PeriodicSquareDistance (1000 times faster) both calculate:
# Delta^2 = \sum_{t,t'} d_{t,t'}^2

######################################################################
# Calculates the square distance between the same configuration at two different times.
# The operation is not vectorized, so it is slow. I used it for debug
# because I tend to mess up with even the simplest vectorized operations.
def ConfSqDistPBC(positions_a, positions_b, box_size):
    sum=0
    for i in range(len(positions_a)):
        sum+=ParticleSqDistPBC(positions_a[i],positions_b[i],box_size)
    return sum

######################################################################
# Calculates the square distance between the same configuration at two different times.
# The operation is vectorized, so use this one instead of ConfDistPBC()
def PeriodicSquareDistance(vec_a, vec_b, box_size):
    delta = np.abs(vec_a - vec_b) #subtraction and abs are done component by component
    delta = np.where(delta > box_size-delta, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return np.square(delta).sum()





######################################################################
# Calculate the mutual distances between all particles

def CalculateRelativeDistances(positions, Natoms, L):
#Takes a vector of n positions.
#Returns a vector of n(n-1)/2 distances
    distances=np.zeros(Natoms*(Natoms-1)/2,dtype=np.double)
    pos=0
    for i in range(Natoms-1):
        distances[pos:pos+Natoms-1-i]+=ParticleDistPBC(positions[i],positions[i+1:],np.array([L,L,L]))
        pos+=Natoms-i-1
    return distances












#-------------------------#
###########################
#      DISPLACEMENTS      #
###########################
#-------------------------#
#-----------------------------------------
# -- Displacement between time t and t' --
#-----------------------------------------
# The functions of this section take the positions of two configurations and return a configuration.

######################################################################
# Periodic displacement 
# Returns a vector of the Natoms 3D displacements (one 3D vector for each particle)
def PeriodicDisplacement(vec_a, vec_b, box_size):
    # First we substract one to the other
    delta = vec_a - vec_b
    #Then we apply periodic boundary conditions through a double ternary
    return np.where(delta > 0.5 * box_size, delta - box_size, np.where(delta < -0.5 * box_size, delta+box_size, delta))

######################################################################
# Periodic intermediate point
# Given two lists of 3D points, it returns a list of points half way in between.
# For two points xA and xB, it is usually (xA+xB)/2, unless they are 
# across the periodic boundary conditions.
#
# Poi devo scrivere questa funzione in modo meno inefficiente
def PeriodicIntermPoints(vec_a, vec_b, L):
    Lm=0.5*L
    delta = np.abs(vec_a - vec_b)
    out = np.where(delta < Lm, 0.5*(vec_a+vec_b), np.where( vec_a+vec_b<0, 0.5*(vec_a+vec_b+L),0.5*(vec_a+vec_b-L)) )
    return out






#-------------------------#
###########################
#        OVERLAPS         #
###########################
#-------------------------#

######################################################################
# Overlap with the snapshots as input
def OverlapConfs(conf1, conf2,box_size):
    posizioni1=np.array(conf1.particles.position, dtype=np.float64)
    posizioni2=np.array(conf2.particles.position, dtype=np.float64)
    return OverlapPos(posizioni1,posizioni2,box_size)

######################################################################
# Overlap with the positions as input
def OverlapPos(posizioni1, posizioni2,box_size):
    disp=PeriodicDisplacement(posizioni1,posizioni2,box_size).sum(axis=1)
    return OverlapDisp(disp,box_size)

######################################################################
# Overlap with the distance between confs as input
def OverlapDisp(disp,box_size):
    delta=0.13333 #1/3 of the smallest particle diameter
    return np.where(disp<delta,1.,0.).sum()/len(disp)




#-------------------------------------#
#######################################
#          STATIC OBSERVABLES         #
#######################################
#-------------------------------------#

######################################################################
# Calculate pair correlation function

def CalculatePairCorrelationFunction(distances, Natoms, dr=0.1,rMax=2, number_density=1.0):
#Calculate pair correlation function from a clean list of the distances
#between particles (distances between (a,b), (b,a) should appear only once)
#Input:
# -A list of Natoms(Natoms-1)/2 positive distances
#Output:
# -values of g(r): gofr
#
    nbins=int(rMax/dr)+1
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














#-------------------------------------#
#######################################
#         DYNAMIC OBSERVABLES         #
#######################################
#-------------------------------------#
 
######################################################################
# Calculate Mean Square Displacement

def CalculateMeanSquareDisplacements(old_pos, new_pos, box_size):
   assert(len(old_pos==65))
   assert(len(old_pos==new_pos))
   return PeriodicSquareDistance(old_pos, new_pos, box_size)/len(old_pos) #sum of the squared displacements divided by Natoms
 

######################################################################
# Self-intermediate scattering function

@jit (nogil=True, nopython=False)
def ComputeFkt(NX, NY, NZ, L, displacements):
    #Given the displacement vector and the wave numbers, we calculate the value of Fkt
   Natoms=len(displacements)
   Fk_Deltat = np.zeros( Natoms, dtype=np.complex128)
   numk=48 #2(+/-) x 3(dimensions) x 6(permutations) = 48
   for nx in [NX,-NX]:
      for ny in [NY,-NY]:
         ## note that we could spare this last loop and double the weight of these z-conributions and take real part.. but it's not elegant.
         for nz in [NZ,-NZ]:
            for k_set_index in range(6):
               k_vector = getKSets_function(nx,ny,nz,k_set_index)
               Fk_Deltat += np.exp( (2.0j*np.pi/L) * np.sum(k_vector*displacements,1) ) 
   ## remember that Fk_Deltat is an array (size Natoms)
   #We should check that the imaginary part of Fk is zero
   return np.real(Fk_Deltat).sum()/(numk*Natoms) 

 


######################################################################
# Permutations of the wave vector indices

@jit (nogil=True, nopython=True)
def getKSets_function(nx,ny,nz, it):
    if it==0:
        return np.array([nx,ny,nz])#,dtype=float)
    if it==1:
        return np.array([nz,nx,ny])#,dtype=float)
    if it==2:
        return np.array([ny,nz,nx])#,dtype=float)
    if it==3:
        return np.array([nx,nz,ny])#,dtype=float)
    if it==4:
        return np.array([nz,ny,nx])#,dtype=float)
    if it==5:
        return np.array([ny,nx,nz])#,dtype=float)


######################################################################
#Calculate the height at which Fk cuts height (1/e is default)
def CalculateTau(Fk, firstframe, lastframe, height=0.36787944117144232159,dt=0.0025,every_forMemory=1):
    print("Calculating TAU from Fk(t)")
    value=Fk[firstframe]
    if value<height:
        print("ERROR: Fk(t) starts at ",Fk[firstframe],", which is lower than height=",height)
        raise SystemExit
    if Fk[lastframe]>height:
        print("ERROR: Fk(t) ends at ",Fk[lastframe],", which is higher than height=",height)
        raise SystemExit
    if lastframe-firstframe<10:
        print("We require at least 10 points to interpolate tau (there are currently ",lastframe-firstframe,").")
        raise SystemExit
    iframe=firstframe
    while(value>height):
        iframe+=1
        value=Fk[iframe]

    #The good value is between iframe-1 and iframe
    #I interpolate cubically in the zone around these two points
    reduced_range_start=iframe-3
    reduced_range_end=iframe+3
    x=range(reduced_range_start,reduced_range_end)
    y=Fk[reduced_range_start:reduced_range_end]
    print("start:",reduced_range_start)
    print("end:",reduced_range_end)
    print("len(Fk):",len(Fk))
    print("len(y):",len(y))
    print("len(x):",len(x))

    interp = interp1d(x, y, kind='cubic')

    def f(x):
        return interp(x)-height

    x0 = brentq(f, iframe-1,iframe)
    tau = x0*every_forMemory*dt

    return tau
