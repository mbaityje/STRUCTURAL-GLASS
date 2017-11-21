#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# calculates its energy through two different methods.
#
# To launch a simulation:
# python T10-Energy.py configuration.gsd
#
# For example:
# python T10-Energy.py sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import kobandersen #Here I have put the Kob-Andersen parameters
import matplotlib.pyplot as plt
from numba import jit

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
# Set potential
#
################################################################
NeighborsListLJ = md.nlist.cell()
print(" *** KApotentialShort *** ")
myLjPair=kobandersen.KApotentialShort(NeighborsListLJ)

################################################################
# 
# Set analyzer
#
################################################################
analyzer_quantities = ['temperature', 'potential_energy', 'kinetic_energy']
analyzer = hoomd.analyze.log(filename='prova.txt', overwrite=True, quantities=analyzer_quantities, period=1)

################################################################
# 
# CALCULATE ENERGY
#
################################################################
groupA = hoomd.group.type(name='a-particles', type='A')
groupB = hoomd.group.type(name='b-particles', type='B')

def PotEn():
    modeT.set_params(dt=1e-10)
    hoomd.run(2)
    U=analyzer.query('potential_energy')
    modeT.set_params(dt=0.0025)
    return U


snapshot=system.take_snapshot(dtype='double')
positions=np.array(snapshot.particles.position) #Now 'positions' is a NatomsxDIM vector,
                                                #storing all the particles' positions

#
# From the integrator
#
modeT=md.integrate.mode_standard(dt=0.000000001)
integratorT=md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)
print("Energy of the system from the analyzer:", PotEn())
print("*****************")

#
#Sum over per-particle energies
#
energy=0
for particle in hoomd.group.all():
    energy+=particle.net_energy
energyA=0
for particle in groupA:
    energyA+=particle.net_energy
energyB=0
for particle in groupB:
    energyB+=particle.net_energy
print("Energy from per-particle energy",energy)
print("EnergyA: ",energyA)
print("EnergyB: ",energyB)
print("*****************")

#
#Memory of the last measurement by the analyzer
#
print("Energy from lj.pair.get_energy: ",myLjPair.get_energy(hoomd.group.all()))
print("EnergyA: ",myLjPair.get_energy(groupA))
print("EnergyB: ",myLjPair.get_energy(groupB))
print("*****************")

#
#Sum over the couplings
#
tags=np.linspace(0,Natoms-1,Natoms, dtype=np.int32)
U=0
for i in range(Natoms):
    U += myLjPair.compute_energy(tags1=np.array(tags[i:i+1]), tags2=np.array(tags[i+1:Natoms]))
print("Energy from lj.pair.compute_energy: ",U)
energyAA=0
for i in range(Natoms):
    if system.particles[i].type=='B':
        continue
    tags1=np.array(tags[i:i+1])
    tags2=np.array(tags[i+1:Natoms])
    del_indices=[]
    for ii in range(len(tags2)):
        if(system.particles[int(tags2[ii])].type=='B'):
            del_indices.append(ii)
    tags2=np.delete(tags2,del_indices)
    energyAA += myLjPair.compute_energy(tags1=tags1, tags2=tags2)
print("EnergyAA: ",energyAA)
energyAB=0
for i in range(Natoms):
    if system.particles[i].type=='B':
        continue
    tags1=np.array(tags[i:i+1])
    tags2=np.array(tags[i+1:Natoms])
    del_indices=[]
    for ii in range(len(tags2)):
        if(system.particles[int(tags2[ii])].type=='A'):
            del_indices.append(ii)
    tags2=np.delete(tags2,del_indices)
    energyAB += myLjPair.compute_energy(tags1=tags1, tags2=tags2)
print("EnergyAB: ",energyAB)
energyBA=0
for i in range(Natoms):
    if system.particles[i].type=='A':
        continue
    tags1=np.array(tags[i:i+1])
    tags2=np.array(tags[i+1:Natoms])
    del_indices=[]
    for ii in range(len(tags2)):
        if(system.particles[int(tags2[ii])].type=='B'):
            del_indices.append(ii)
    tags2=np.delete(tags2,del_indices)
    energyBA += myLjPair.compute_energy(tags1=tags1, tags2=tags2)
print("EnergyBA: ",energyBA)
energyBB=0
for i in range(Natoms):
    if system.particles[i].type=='A':
        continue
    tags1=np.array(tags[i:i+1])
    tags2=np.array(tags[i+1:Natoms])
    del_indices=[]
    for ii in range(len(tags2)):
        if(system.particles[int(tags2[ii])].type=='A'):
            del_indices.append(ii)
    tags2=np.delete(tags2,del_indices)
    energyBB += myLjPair.compute_energy(tags1=tags1, tags2=tags2)
print("EnergyBB: ",energyBB)
print("EnergyAA+EnergyAB+EnergyBA+EnergyBB = ",energyAA+energyAB+energyBA+energyBB)
print("*****************")



print("REMEMBER, NONE OF THESE METHODS WORKS IF THE ANALYZER IS NOT SET (MAKING THEM DE FACTO USELESS IN SOME SITUATIONS).")
print("ALSO, THE INFORMATION THEY YIELD DEPENDS ON THE STATE OF THE INTEGRATOR")






print("\nNow calculate the energy from the positions.")
print("These functions are easily optimized, in case one actually needed them")
#XPLOR smoothing function
def S(r, r_on, r_cut):
    assert(r>0)
    if r<r_on:
        return 1
    elif r>r_cut:
        return 0
    else:
        rc2=r_cut*r_cut
        ro2=r_on*r_on
        r2=r*r
        term1=rc2 - r2
        term2=rc2 + 2*r2 - 3*ro2
        term3=rc2 - ro2
        return (term1*term1*term2)/(term3*term3*term3)
#LJ pair pure potential
def Vpair(r, eps, sigma, r_cut):
    assert(r>0)
    if r>=r_cut:
        return 0
    else:
        sigma_on_r=sigma/r
        temp=sigma_on_r*sigma_on_r
        sigma_on_r6=temp*temp*temp
        return 4*eps*(sigma_on_r6*sigma_on_r6 - sigma_on_r6)
#LJ potential with the XPLOR smoothing
def Vij(r, eps, sigma, r_on, r_cut):
    if(r_on<r_cut):
        return S(r,r_on,r_cut)*Vpair(r,eps, sigma,r_cut)
    elif(r_on>=r_cut):
        return Vpair(r,eps,sigma,r_cut)-Vpair(r_cut,eps,sigma,r_cut)
#Parametri del potenziale
#AA
assert(myLjPair.pair_coeff.get_metadata()[1]['typei']=='A')
assert(myLjPair.pair_coeff.get_metadata()[1]['typej']=='A')
eps_AA=myLjPair.pair_coeff.get_metadata()[1]['epsilon']
sig_AA=myLjPair.pair_coeff.get_metadata()[1]['sigma']
rcut_AA=myLjPair.pair_coeff.get_metadata()[1]['r_cut']
ron_AA=myLjPair.pair_coeff.get_metadata()[1]['r_on']
#AB
assert(myLjPair.pair_coeff.get_metadata()[0]['typei']=='A')
assert(myLjPair.pair_coeff.get_metadata()[0]['typej']=='B')
eps_AB=myLjPair.pair_coeff.get_metadata()[0]['epsilon']
sig_AB=myLjPair.pair_coeff.get_metadata()[0]['sigma']
rcut_AB=myLjPair.pair_coeff.get_metadata()[0]['r_cut']
ron_AB=myLjPair.pair_coeff.get_metadata()[0]['r_on']
#BB
assert(myLjPair.pair_coeff.get_metadata()[2]['typei']=='B')
assert(myLjPair.pair_coeff.get_metadata()[2]['typej']=='B')
eps_BB=myLjPair.pair_coeff.get_metadata()[2]['epsilon']
sig_BB=myLjPair.pair_coeff.get_metadata()[2]['sigma']
rcut_BB=myLjPair.pair_coeff.get_metadata()[2]['r_cut']
ron_BB=myLjPair.pair_coeff.get_metadata()[2]['r_on']
#Calcolo dell'energia per gruppi
energiaAA=0
for p1 in groupA:
    for p2 in groupA:
        if p1.tag>=p2.tag:
            continue
        r=PeriodicDistance(np.array(p1.position),np.array(p2.position),L)
        energiaAA+=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)
energiaAB=0
for p1 in groupA:
    for p2 in groupB:
        r=PeriodicDistance(np.array(p1.position),np.array(p2.position),L)
        energiaAB+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
energiaBB=0
for p1 in groupB:
    for p2 in groupB:
        if p1.tag>=p2.tag:
            continue
        r=PeriodicDistance(np.array(p1.position),np.array(p2.position),L)
        energiaBB+=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)

print("energiaAA= ", energiaAA)
print("energiaAB= ", energiaAB)
print("energiaBB= ", energiaBB)
print("Energy measured manually from the positions: ", energiaAA+energiaBB+energiaAB)
print("*****************")





print("\nNow calculate the energy from a snapshot")
#Calcolo dell'energia per gruppi
energiaAA=0
for p1 in groupA:
    for p2 in groupA:
        if p1.tag>=p2.tag:
            continue
        pos1=np.array(snapshot.particles.position[p1.tag])
        pos2=np.array(snapshot.particles.position[p2.tag])
        r=PeriodicDistance(pos1,pos2,L)
        energiaAA+=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)
energiaAB=0
for p1 in groupA:
    for p2 in groupB:
        pos1=np.array(snapshot.particles.position[p1.tag])
        pos2=np.array(snapshot.particles.position[p2.tag])
        r=PeriodicDistance(pos1,pos2,L)
        energiaAB+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
energiaBB=0
for p1 in groupB:
    for p2 in groupB:
        if p1.tag>=p2.tag:
            continue
        pos1=np.array(snapshot.particles.position[p1.tag])
        pos2=np.array(snapshot.particles.position[p2.tag])
        r=PeriodicDistance(pos1,pos2,L)
        energiaBB+=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)
print("energiaAA= ", energiaAA)
print("energiaAB= ", energiaAB)
print("energiaBB= ", energiaBB)
print("Energy measured manually from the positions: ", energiaAA+energiaBB+energiaAB)
print("*****************")
