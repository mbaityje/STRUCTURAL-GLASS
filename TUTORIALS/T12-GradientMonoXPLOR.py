#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################
#
# DESCRIPTION This example reads a configuration from a .gsd file, and
# calculates the gradient of the energy with a Kob-Andersen
# monodisperse potential smoothened through XPLOR.
#
# Then the energy is minimized, and the gradient is calculated again,
# just to show that it is now zero.  Finally, the initial
# configuration is taken again, and instead of minimizing the energy
# we minimize the square of the gradient.
#
# To launch a simulation:
# python T12-GradientMonoXPLOR.py configuration.gsd
#
# For example:
# python T12-GradientMonoXPLOR.py sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd
# python T12-GradientMonoXPLOR.py sample-states/N65.gsd
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

import time
def profile(func, fargs):
    start = time.time()
    output=func(*fargs)
    end = time.time()
    print("Time elapsed: ",end - start," s")
    return output

def PeriodicDistance(vec_a, vec_b, box_size):
# This function measures the distances between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #substraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return np.sqrt((delta ** 2).sum(axis=-1))

def PeriodicDisplacement(vec_a, vec_b, box_size):
    # First we substract one to the other
    delta = vec_a - vec_b
    #Then we apply periodic boundary conditions through a double ternary
    return np.where(delta > 0.5 * box_size, delta - box_size, np.where(delta < -0.5 * box_size, delta+box_size, delta))

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
eps=1
sigma=1
r_on=1.2
r_cut=1.4
print(" *** Creating a monodisperse LJ system *** ")
NeighborsListLJ = md.nlist.cell()
NeighborsListLJ.set_params(r_buff=0.0)
myLjPair = md.pair.lj(r_cut=r_cut, nlist=NeighborsListLJ)
myLjPair.pair_coeff.set('A', 'A', epsilon=eps, sigma=sigma, r_cut=r_cut, r_on=r_on)
myLjPair.pair_coeff.set('A', 'B', epsilon=eps, sigma=sigma, r_cut=r_cut, r_on=r_on)
myLjPair.pair_coeff.set('B', 'B', epsilon=eps, sigma=sigma, r_cut=r_cut, r_on=r_on)
myLjPair.set_params(mode="xplor")

################################################################
# 
# Set analyzer
#
################################################################
analyzer_quantities = ['temperature', 'potential_energy', 'kinetic_energy']
analyzer = hoomd.analyze.log(filename=None, overwrite=True, quantities=analyzer_quantities, period=1)


################################################################
# 
# DEFINE SOME MORE FUNCTIONS REGARDING THIS TUTORIAL
#
################################################################

#Returns system's energy if the analyzer is set
def PotEn():
    modeT.set_params(dt=1e-12)
    hoomd.run(2)
    U=analyzer.query('potential_energy')
    modeT.set_params(dt=0.0025)
    return U

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
    
#Derivative of the LJ pair pure potential
#I include an extra 1/r factor that would need to be put later on
#because of the chain rule derivation dr/d\vec{r}
def VpairPrime(r, eps, sigma, r_cut):
    if r>=r_cut:
        return 0;
    invr=1./r
    invr2=invr*invr
    invr4=invr2*invr2
    invr6=invr4*invr2
    invr8=invr4*invr4
    s2=sigma*sigma
    s6=s2*s2*s2
    return 48*eps*s6*invr8*(0.5 - s6*invr6)

#Calculate energy explicitly from the snapshot.
#Instead of calculating the energy of each group, I go through all the particles,
#and only later I check which group they belong to.
def CalculateSnapEnergySlower():
    energia=0
    for p1 in hoomd.group.all():
        for p2 in hoomd.group.all():
            if p1.tag>=p2.tag:
                continue
            pos1=np.array(snapshot.particles.position[p1.tag], dtype=np.float64)
            pos2=np.array(snapshot.particles.position[p2.tag], dtype=np.float64)
            r=PeriodicDistance(pos1,pos2,L)
            energia+=Vpair(r,eps,sigma,r_cut)
    return energia
def S(r, r_on, r_cut):
    assert(r>0)
    if r<r_on:
        return 1
    elif r>r_cut:
        return 0
    else:
        return Stilde(r, r_on, r_cut)
#Stilde is the non-trivial part of S
def Stilde(r, r_on, r_cut):
    assert(r>r_on)
    assert(r<r_cut)
    rc2=r_cut*r_cut
    ro2=r_on*r_on
    r2=r*r
    term1=rc2 - r2
    term2=rc2 + 2*r2 - 3*ro2
    term3=rc2 - ro2
    return (term1*term1*term2)/(term3*term3*term3)

#For the derivative of S, I only consider the non-trivial part
def StildePrime(r, r_on, r_cut):
    assert(r>r_on)
    assert(r<r_cut)
    rc2=r_cut*r_cut
    ro2=r_on*r_on
    r2=r*r
    term1=rc2 - r2
    term2=ro2 - r2
    term3=rc2 - ro2
    return 12*r*term1*term2/(term3*term3*term3)

#LJ potential with the XPLOR smoothing
def Vij(r, eps, sigma, r_on, r_cut):
    if r>=r_cut:
        return 0
    if(r_on<r_cut):
        return S(r,r_on,r_cut)*Vpair(r,eps, sigma,r_cut)
    elif(r_on>=r_cut):
        return Vpair(r,eps,sigma,r_cut)-Vpair(r_cut,eps,sigma,r_cut)

#Derivative of the potential with XPLOR smoothing
def Vprime(r, eps, sigma, r_on, r_cut):
    if r>=r_cut:
        return 0
    elif r<=r_on:
        return VpairPrime(r, eps, sigma, r_cut)
    else:
        return VpairPrime(r, eps, sigma, r_cut)*Stilde(r,r_on,r_cut)+(Vpair(r, eps, sigma, r_cut)*StildePrime(r, r_on, r_cut))/r

#Minimize the energy of a snapshot
def Minimize(snap):
    system.restore_snapshot(snap)
    fire.cpp_integrator.reset()
    while not(fire.has_converged()):
        hoomd.run(100)       
    eIS=analyzer.query('potential_energy')
    return eIS

#Calculate Gradient of the LJ potential
def CalculateGradient(snapshot):
    energia=0
    gradVector = np.zeros((Natoms, 3))
    for i in range(Natoms):
        for j in range(i+1,Natoms):
            posi=np.array(snapshot.particles.position[i], dtype=np.float64)
            posj=np.array(snapshot.particles.position[j], dtype=np.float64)
            rvec=PeriodicDisplacement(posi,posj,L)
            r=np.sqrt(np.square(rvec).sum())

            energia+=Vij(r, eps, sigma, r_on, r_cut)
            Vp=Vprime(r, eps, sigma, r_on, r_cut)
            temp=rvec*Vp

            gradVector[i] +=  temp 
            gradVector[j] += -temp
    return energia/Natoms,gradVector/Natoms


################################################################
# 
# CALCULATE ENERGY AND GRADIENT
#
################################################################

modeT=md.integrate.mode_standard(dt=1e-12)
integratorT=md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)


# Energy from the integrator
eAnalyzer=profile(PotEn,())/Natoms
integratorT.disable()
print("e:", eAnalyzer)
F=np.array([system.particles[i].net_force for i in range(Natoms)])/Natoms
F2=np.square(F).sum()
print("F2= ",F2)

snapshot=system.take_snapshot(dtype='double')
energia,gradVector=CalculateGradient(snapshot)
G2=np.square(gradVector).sum()
print("e=",energia)
print("G2=",G2)




#Minimize the energy
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.5, ftol=1e-3, Etol=1e-8, wtol=1e-3)
integrator = md.integrate.nve(group=hoomd.group.all())
while not(fire.has_converged()):
   hoomd.run(100)
eIS=Minimize(snapshot)/Natoms
print("\nMinimize Energy...\neIS=",eIS)

#Calculate gradient again
snapnew=system.take_snapshot(dtype='double')
F=np.array([system.particles[i].net_force for i in range(Natoms)])/Natoms
F2=np.square(F).sum()
print("F2= ",F2)

snapshot=system.take_snapshot(dtype='double')
energia,gradVector=CalculateGradient(snapshot)
G2=np.square(gradVector).sum()
print("e=",energia)
print("G2=",G2)


raise SystemExit

#NOTE:
#Per calcolare la forza dall'integrator, posso sia prenderla dal LJpair, che dalla feature della particella
#Per esempio i dati delle forze di LJ sulla particella 0 si ottengono via:
myLJpair.forces[0] #Non so come si spacchetta per ottenere solo la forza
p = system.particles[0]
p.net_force

#Nello snapshot non vedo la possibilita` di prendere le forze, ma se
#le particelle hanno massa unitaria (come nel mio caso), allora posso usare
acc=snapshot.particles.acceleration

