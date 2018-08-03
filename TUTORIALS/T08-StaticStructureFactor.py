#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# calculates the static structure factor S(k).
# In addition, the pair correlation function is also calculated,
# in order to perform a check.
#
# A VERY NICE REVIEW ON S(k):    arXiv:1606.03610
#
# To launch a simulation:
# python T8-StaticStructureFactor.py configuration.gsd
#
# For example:
# python T8-StaticStructureFactor.py sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import kobandersen #Here I have put the Kob-Andersen parameters
import matplotlib.pyplot as plt
from numba import jit

Nmax=30 #Maximum wave vector: (2pi/L)[Nmax,Nmax,Nmax]

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

def CalculateRelativeDistanceVectors(positions):
#Takes a vector of n positions.
#Returns a vector of n(n-1)/2 distance vectors
# NO PERIODIC BOUNDARY CONDITIONS ARE APPLIED
    rij=np.ndarray(shape=(Natoms*(Natoms-1)/2,3),dtype=np.double)
    pos=0
    for i in range(Natoms-1):
        rij[pos:pos+Natoms-1-i]+=(positions[i]-positions[i+1:])
        pos+=Natoms-i-1
    return rij


def CalculateAllRelativeDistances(positions):
#Takes a vector of n positions.
#Returns a vector of n(n-1)/2 distances
    distances=np.zeros(Natoms*Natoms,dtype=np.double)
    for i in range(Natoms-1):
        for j in range(Natoms-1):
            distances[i+Natoms*j]+=PeriodicDistance(positions[i],positions[j],np.array([Lx,Ly,Lz]))
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

@jit (nogil=True, nopython=False)
def CalculateStaticStructureFactorFromPositions(positions, L, dk=-1, NXmax=-1, NYmax=-1, NZmax=-1, Nmax=3):
#
# S(k) = <[sum_i cos(k·r_i)]^2+[sum_i sin(k·r_i)]^2>
#
#Calculate the Static Structure Factor from a clean list of all the positions
#Input:
# -list of Natoms positions
# -linear system size
# -delta k bins
# -maximum wave numbers for measurements
#Output:
# -values of S(k): Sk
#
    if NXmax<0:
        NXmax=Nmax
    if NYmax<0:
        NYmax=Nmax
    if NZmax:
        NZmax=Nmax
    if dk<0:
        dk=2*np.pi/L

    factor=2*np.pi/L
    kMin=factor
    kMax=factor*np.linalg.norm([NXmax,NYmax,NZmax])
    assert(kMin>=0 and kMax>kMin) #Usually kMin=2pi/L, kMax=2pi/(smallest particle diameter)
    assert(dk>0)
    nbins=int((kMax-kMin)/dk)+1

    Sk = np.zeros(nbins)
    ninbin = np.zeros(nbins)
    klist = np.ndarray(nbins)
    
    for NX in range(NXmax):
        for NY in range(NYmax):
            for NZ in range(NZmax):
                value=ComputeSk(NX, NY, NZ, L, positions)
                k=factor*np.linalg.norm([NX,NY,NZ])
                ibin=int((k-kMin)/dk)
                if (ibin>=0 and ibin<nbins):
                    Sk[ibin]+=value
                    ninbin[ibin]+=1
    for ibin in range(nbins):
        klist[ibin]=kMin+(ibin+0.5)*dk
        if ninbin[ibin]>0:
            Sk[ibin]/=ninbin[ibin]
        else:
            Sk[ibin]=np.nan
    return klist,Sk

@jit (nogil=True, nopython=False)
def ComputeSk(NX, NY, NZ, L, positions):
    Natoms=np.double(len(positions))
    factor=2*np.pi/L
    output=0
    for nx in [NX,-NX]:
        for ny in [NY,-NY]:
            for nz in [NZ,-NZ]:
                for k_set_index in range(6):
                    n_vector = getKSets_function(nx,ny,nz,k_set_index)
                    k=factor*n_vector
                    kr=np.inner(k,positions) #kr is an Natoms-dimensional vector
                    coseno=np.cos(kr)
                    seno=np.sin(kr)
                    output+=(coseno.sum()**2+seno.sum()**2)
    return np.double(output)/(48*Natoms)

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
# CALCULATE STATIC STRUCTURE FACTOR
#
################################################################

snapshot=system.take_snapshot(dtype='double')
positions=np.array(snapshot.particles.position) #Now 'positions' is a NatomsxDIM vector,
                                      #storing all the particles' positions

#First Method, from positions (has a bias when dealing with a single sample)
klist,Sk=CalculateStaticStructureFactorFromPositions(positions, L, dk=2*np.pi/L, Nmax=Nmax)

#Second Method, from distances
#distances=CalculateRelativeDistanceVectors(positions)
#klist,Sk=CalculateStaticStructureFactorFromDistances(distances, L, Natoms, dk=2*np.pi/L, Nmax=Nmax)

################################################################
# 
# CALCULATE PAIR CORRELATION FUNCTION
#
################################################################

#ACCESS PARTICLE POSITIONS
#To access the array of the particle positions I need to take a
#snapshot (don't ask why it can't be accessed from the system
#directly, but it can't). Anyhow, it's just one more line of code.

#I set a couple of parameters for the binning of the g(r)
rMax=L/2
dr=0.01
distances=CalculateRelativeDistances(positions)
gofr=CalculatePairCorrelationFunction(distances, Natoms, dr=dr, rMax=rMax, number_density=rho)
rvalues = np.arange(0,rMax,dr)
r=rvalues-0.5*dr


#As a check, we calculate the Static Structure Factor by antitransforming the g(r)
nbins=np.int(len(klist))
factor=np.double(4*np.pi*rho)
Sk_check=np.ndarray(nbins)
hofr=np.array(gofr-1)

for i in range(0,nbins):
    k=klist[i]
    kr=np.array(k*r)
    somma=np.sum(factor*dr*hofr*(r**2)*np.sin(kr)/kr)
    Sk_check[i]=1+somma
    
################################################################
# 
# FIGURES
# 
################################################################
#S(k)
output_Sk=np.column_stack((klist, Sk))
output_Skcheck=np.column_stack((klist, Sk_check))
np.savetxt('./test-output/Sk.txt',output_Sk,fmt='%g %.14g')
np.savetxt('./test-output/Skcheck.txt',output_Skcheck,fmt='%g %.14g')
plt.plot(klist, Sk, label='S(k)')
plt.plot(klist, Sk_check, label='S(k) from g(r)')
plt.xlabel('k')
plt.ylabel('S(k)')
plt.title('Static Structure Factor of '+filename)
plt.grid(True)
plt.savefig(filename+"_Sk.png")
plt.show()

#g(r)
output_gofr=np.column_stack((rvalues, gofr))
np.savetxt('./test-output/gofr.txt',output_gofr,fmt='%g %.14g')
plt.plot(rvalues, gofr)
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Pair Correlation function of '+filename)
plt.grid(True)
plt.savefig(filename+"_PAIR.png")
plt.show()











