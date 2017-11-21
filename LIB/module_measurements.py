#!/usr/bin/python
import numpy as np
from numba import jit
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy import optimize as opt



######################################################################
# Function to check whether the sample reached equilibrium

def IsAtEquilibrium34(array, nSigma=1, rtol=0.01, maxErr=0.005, maxStd=0.1, verbose=False):
#Function takes observable in list, and compares its value in the
# third and fourth quarter of the list. Returns True if the two are compatible.
# Conditions for equilibration:
# -The two averages are statistically compatible
# -The two averages are within 1% relative difference
# -The statistical error of the fourth quarter is relatively <1%
# -The standard deviation of the fourth quarter is relatively <10% 
   quarter4=len(array)
   quarter3=(int)(3.*quarter4/4)
   quarter2=(int)(quarter4/2)
   part3=array[quarter2:quarter3]
   part4=array[quarter3:quarter4]
   mean3=part3.mean()
   mean4=part4.mean()
   std4=part4.std()
   err3=part3.std()/np.sqrt(quarter3-quarter2-1) #Statistical error of third quarter
   err4=std4/np.sqrt(quarter4-quarter3-1) #Statistical error of fourth quarter
   compatibility= abs(mean3-mean4)/np.sqrt(err3*err3+err4*err4) #Compatibility of the two measurements (in units of errorbars)
   relative_distance= abs((mean3-mean4)/mean4)
   compatible= True if compatibility<nSigma else False #nSigma std err apart
   relatively_close = True if relative_distance<rtol else False #1% tolerance
   small_errorbar = True if err4/abs(mean4)<maxErr else False
   small_std = True if std4/abs(mean4)<maxStd else False
   if(verbose):
      print("quarter4: ",quarter4,"quarter3: ",quarter3,"quarter2: ",quarter2)
      print("Mean3: ", mean3," Err3: ",err3)
      print("Mean4: ", mean4," Err4: ",err4)
      print("Compatibility: ",compatibility," --> Compatible=",compatible)
      print("Relative Distance: ",relative_distance," --> Relatively close=",relatively_close)
      print("Relative Error Bar: ", err4/abs(mean4)," --> Small Error Bar=",small_errorbar)
      print("Relative Standard Deviation: ", std4/abs(mean4)," --> Small Standard Dev.=",small_std)
   return (compatible and relatively_close and small_errorbar and small_std) #All checks must succeed



######################################################################
# Calculate periodic square distance between two points

def PeriodicSquareDistance(vec_a, vec_b, box_size):
# This function measures the distance between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #subtraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return (delta ** 2).sum(axis=-1)

######################################################################
# Calculate periodic distance between two points

def PeriodicDistance(vec_a, vec_b, box_size):
    return np.sqrt(PeriodicSquareDistance(vec_a, vec_b, box_size))

######################################################################
# Calculate periodic distance between two points, in another way

def PeriodicDifference(vec_a, vec_b, box_size):
# This function measures the distance between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #substraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return delta.sum(axis=-1)

######################################################################
# Calculate Mean Square Displacement

def CalculateMeanSquareDisplacements(old_pos, new_pos, box_size):
   all_displacements=PeriodicSquareDistance(old_pos, new_pos, box_size)
   return all_displacements.sum()/len(all_displacements) #sum of the squared displacements divided by Natoms
 
######################################################################
# Calculate the mutual distances between all particles

def CalculateRelativeDistances(positions, Natoms, L):
#Takes a vector of n positions.
#Returns a vector of n(n-1)/2 distances
    distances=np.zeros(Natoms*(Natoms-1)/2,dtype=np.double)
    pos=0
    for i in range(Natoms-1):
        distances[pos:pos+Natoms-1-i]+=PeriodicDistance(positions[i],positions[i+1:],np.array([L,L,L]))
        pos+=Natoms-i-1
    return distances

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


 
######################################################################
# Self-intermediate scattering function

@jit (nogil=True, nopython=False)
def ComputeFkt(NX, NY, NZ, L, displacements):
    #Given the displacement vector and the wave numbers, we calculate the value of Fkt
   Natoms=len(displacements)
   Fk_Deltat = np.zeros( len(displacements), dtype=np.complex128)
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
# Periodic displacement (maintains information on the sign)

def PeriodicDisplacement(vec_a, vec_b, box_size):
# This function measures the vector displacement between two lists of positions, vec_a and vec_b,
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    # First we substract one to the other
    delta = vec_a - vec_b
    #Then we apply periodic boundary conditions through a double ternary
    return np.where(delta > 0.5 * box_size, delta - box_size, np.where(delta < -0.5 * box_size, delta+box_size, delta))


######################################################################
# Overlap with the configurations as input
def OverlapConfs(conf1, conf2,box_size):
   posizioni1=np.array(conf1.particles.position)
   posizioni2=np.array(conf2.particles.position)
   return OverlapPos(posizioni1,posizioni2,box_size)

######################################################################
# Overlap with the positions as input
def OverlapPos(posizioni1, posizioni2,box_size):
   dist=PeriodicDifference(posizioni1,posizioni2,box_size)
   return OverlapDist(dist,box_size)

######################################################################
# Overlap with the distance between confs as input
def OverlapDist(dist,box_size):
   delta=0.2 #less than half small particle diameter
   return np.where(dist<delta,1.,0.).sum()/len(dist)

######################################################################
#Calculate the height at which Fk cuts height (1/e is default)
def CalculateTau(Fk, firstframe, lastframe, height=0.36787944117144232159,dt=0.0025,every_forMemory=1):
    print("Calculating TAU from Fk(t)")
    value=Fk[firstframe]
    if value<height:
        print("ERROR: Fk(t) starts from lower than height=",height)
        sys.exit()
    if Fk[lastframe]>height:
        print("ERROR: Fk(t) ends at higher than height=",height)
        sys.exit()
    if lastframe-firstframe<10:
        print("We require at least 10 points to interpolate tau.")
        sys.exit()
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
    interp = interp1d(x, y, kind='cubic')

    def f(x):
        return interp(x)-height

    x0 = brentq(f, iframe-1,iframe)
    tau = x0*every_forMemory*dt

    return tau
