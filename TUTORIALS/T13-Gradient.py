#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################
#
# DESCRIPTION This example reads a configuration from a .gsd file, and
# calculates the gradient of the energy with a bidisperse Kob-Andersen
# potential.  Then the energy is minimized, and the gradient is
# calculated again, just to show that it is now zero.  Finally, the
# initial configuration is taken again, and instead of minimizing the
# energy we minimize the square of the gradient.
#
# We compare the gradient gradU with the vector of the forces F in the
# integrator.  We should have F=-gradU.
#
# To launch a simulation:
# python T13-Gradient.py configuration.gsd
#
# For example:
# python T13-Gradient.py sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd
# python T13-Gradient.py sample-states/N65.gsd
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
NeighborsListLJ = md.nlist.cell()
print(" *** KApotentialShort *** ")
myLJpair=kobandersen.KApotentialShort(NeighborsListLJ)

#Parametri del potenziale
#AA
assert(myLJpair.pair_coeff.get_metadata()[1]['typei']=='A')
assert(myLJpair.pair_coeff.get_metadata()[1]['typej']=='A')
eps_AA=myLJpair.pair_coeff.get_metadata()[1]['epsilon']
sig_AA=myLJpair.pair_coeff.get_metadata()[1]['sigma']
rcut_AA=myLJpair.pair_coeff.get_metadata()[1]['r_cut']
ron_AA=myLJpair.pair_coeff.get_metadata()[1]['r_on']
#AB
assert(myLJpair.pair_coeff.get_metadata()[0]['typei']=='A')
assert(myLJpair.pair_coeff.get_metadata()[0]['typej']=='B')
eps_AB=myLJpair.pair_coeff.get_metadata()[0]['epsilon']
sig_AB=myLJpair.pair_coeff.get_metadata()[0]['sigma']
rcut_AB=myLJpair.pair_coeff.get_metadata()[0]['r_cut']
ron_AB=myLJpair.pair_coeff.get_metadata()[0]['r_on']
#BB
assert(myLJpair.pair_coeff.get_metadata()[2]['typei']=='B')
assert(myLJpair.pair_coeff.get_metadata()[2]['typej']=='B')
eps_BB=myLJpair.pair_coeff.get_metadata()[2]['epsilon']
sig_BB=myLJpair.pair_coeff.get_metadata()[2]['sigma']
rcut_BB=myLJpair.pair_coeff.get_metadata()[2]['r_cut']
ron_BB=myLJpair.pair_coeff.get_metadata()[2]['r_on']

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
    modeT.set_params(dt=1e-13)
    hoomd.run(2)
    U=analyzer.query('potential_energy')
    modeT.set_params(dt=0.0025)
    return U

#XPLOR smoothing function
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
def VpairPrime(r, eps, sigma, r_cut):
    if r>=r_cut:
        return 0
    invr=1./r
    invr2=invr*invr
    invr4=invr2*invr2
    invr6=invr4*invr2
    invr8=invr4*invr4
    s2=sigma*sigma
    s6=s2*s2*s2
    return 48*eps*s6*invr8*(0.5 - s6*invr6)


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

#Calculate energy explicitly from the positions
def CalculatePosEnergySlow():
    energiaAA=0
    for p1 in groupA:
	for p2 in groupA:
	    if p1.tag>=p2.tag:
	        continue
	    r=PeriodicDistance(np.array(p1.position, dtype=np.float64),np.array(p2.position, dtype=np.float64),L)
	    energiaAA+=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)

    energiaAB=0
    for p1 in groupA:
	for p2 in groupB:
	    r=PeriodicDistance(np.array(p1.position, dtype=np.float64),np.array(p2.position, dtype=np.float64),L)
	    energiaAB+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
    energiaBB=0
    for p1 in groupB:
	for p2 in groupB:
	    if p1.tag>=p2.tag:
	        continue
	    r=PeriodicDistance(np.array(p1.position, dtype=np.float64),np.array(p2.position, dtype=np.float64),L)
	    energiaBB+=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)
    return (energiaAA+energiaBB+energiaAB)/Natoms

#Calculate energy explicitly from the snapshot
def CalculateSnapEnergySlow(snapshot):
    energiaAA=0
    for p1 in groupA:
        for p2 in groupA:
            if p1.tag>=p2.tag:
                continue
            pos1=np.array(snapshot.particles.position[p1.tag], dtype=np.float64)
            pos2=np.array(snapshot.particles.position[p2.tag], dtype=np.float64)
            r=PeriodicDistance(pos1,pos2,L)
            energiaAA+=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)
    energiaAB=0
    for p1 in groupA:
        for p2 in groupB:
            pos1=np.array(snapshot.particles.position[p1.tag], dtype=np.float64)
            pos2=np.array(snapshot.particles.position[p2.tag], dtype=np.float64)
            r=PeriodicDistance(pos1,pos2,L)
            energiaAB+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)

    energiaBB=0
    for p1 in groupB:
        for p2 in groupB:
            if p1.tag>=p2.tag:
                continue
            pos1=np.array(snapshot.particles.position[p1.tag], dtype=np.float64)
            pos2=np.array(snapshot.particles.position[p2.tag], dtype=np.float64)
            r=PeriodicDistance(pos1,pos2,L)
            energiaBB+=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)
    return (energiaAA+energiaBB+energiaAB)/Natoms

#Calculate energy explicitly from the snapshot.
#Instead of calculating the energy of each group, I go through all the particles,
#and only later I check which group they belong to.
def CalculateSnapEnergySlower(snapshot):
    energia=0
    for p1 in hoomd.group.all():
        for p2 in hoomd.group.all():
            if p1.tag>=p2.tag:
                continue
            pos1=np.array(snapshot.particles.position[p1.tag], dtype=np.float64)
            pos2=np.array(snapshot.particles.position[p2.tag], dtype=np.float64)
            r=PeriodicDistance(pos1,pos2,L)
            if p1.type=='A':
                if p2.type=='A':
                    energia+=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)
                elif p2.type=='B':
                    energia+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
                else:
                    print("Particle type should be only A or B.")
                    raise SystemExit
            elif p1.type=='B':
                if p2.type=='A':
                    energia+=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
                elif p2.type=='B':
                    energia+=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)
            else:
                print("Particle type should be only A or B.")
                raise SystemExit
    return energia/Natoms

#Minimize the energy of a snapshot
def Minimize(snap):
    system.restore_snapshot(snap)
    fire.cpp_integrator.reset()
    while not(fire.has_converged()):
        hoomd.run(100)       
    eIS=analyzer.query('potential_energy')
    return eIS/Natoms

#Calculate Gradient of the smoothened LJ potential
def CalculateGradient(snapshot):
    energia=0
    gradVector = np.zeros((Natoms, 3))
    for i in range(Natoms):
        for j in range(i+1,Natoms):
            posi=np.array(snapshot.particles.position[i], dtype=np.float64)
            posj=np.array(snapshot.particles.position[j], dtype=np.float64)
            typei=snapshot.particles.typeid[i]
            typej=snapshot.particles.typeid[j]
            rvec=PeriodicDisplacement(posi,posj,L)
            r=np.sqrt(np.square(rvec).sum())

            #See what group the particles belong to
            if typei==typej:
                if typej==0:
                    if r>rcut_AA:
                        continue
                    V=Vij(r,eps_AA,sig_AA,ron_AA,rcut_AA)
                    Vp=Vprime(r, eps_AA, sig_AA, ron_AA, rcut_AA)
                else:
                    if r>rcut_BB:
                        continue
                    Vp=Vprime(r, eps_BB, sig_BB, ron_BB, rcut_BB)
                    V=Vij(r,eps_BB,sig_BB,ron_BB,rcut_BB)
            else:
                if r>rcut_AB:
                    continue
                Vp=Vprime(r, eps_AB, sig_AB, ron_AB, rcut_AB)
                V=Vij(r,eps_AB,sig_AB,ron_AB,rcut_AB)
            temp=rvec*Vp
            energia+=V
            gradVector[i] +=  temp 
            gradVector[j] += -temp
                
    return energia/Natoms,gradVector/Natoms


################################################################
# 
# CALCULATE ENERGY
#
################################################################

groupA = hoomd.group.type(name='a-particles', type='A')
groupB = hoomd.group.type(name='b-particles', type='B')
groupAll=hoomd.group.all()

snapshot=system.take_snapshot(dtype='double')
positions=np.array(snapshot.particles.position, dtype=np.float64) #Now 'positions' is a NatomsxDIM vector,
                                                #storing all the particles' positions
modeT=md.integrate.mode_standard(dt=1e-13)
integratorT=md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)


#
# First we calculate the energy manually and through the integrator,
# to make sure that parameters of the potential we use are the same
# that are implemented in the analyzer.
# Hence, the first part of this tutorial comes directly from tutorial T10.
#






# Energy from the integrator
eAnalyzer=PotEn()/Natoms
integratorT.disable()
print("e:", eAnalyzer)
#Calculate gradient from forces
F=np.array([system.particles[i].net_force for i in range(Natoms)])/Natoms
F2=np.square(F).sum()
print("F2= ",F2)

#Calculate gradient manually
e,G=CalculateGradient(snapshot)
G2=np.square(G).sum()
print("e:", e)
print("G2= ",G2)


#Minimize the energy
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.99, ftol=1e-5, Etol=1e-12, wtol=1e-5)
integrator = md.integrate.nve(group=hoomd.group.all())
while not(fire.has_converged()):
   hoomd.run(100)
eIS=Minimize(snapshot)
print("\n\neIS=",eIS)

#Calculate gradient again
F=np.array([system.particles[i].net_force for i in range(Natoms)])/Natoms
F2=np.square(F).sum()
print("F2= ",F2)


snapnew=system.take_snapshot(dtype='double')
eIS,G=CalculateGradient(snapnew)
g2new=np.square(G).sum()
print("eIS=",eIS)
print("G2=",g2new)

#eISbis=CalculateSnapEnergySlow(snapnew)
#print("eISbis=",eISbis)

raise SystemExit

############################################################
#                                                          #
# Just for fun, let us draw some features of the potential #
#                                                          #
############################################################
import matplotlib.pyplot as plt
#
# XPLOR smoothin function S(r) and its derivative S'(r)
#
plt.figure(1)
plt.xlabel('$r$')
plt.ylabel('$S(r) or S\'(r)$')
#AA
t = np.arange(ron_AA+0.0001, rcut_AA-0.0001, (rcut_AA-ron_AA)/500)
St=[Stilde(i,ron_AA,rcut_AA) for i in t]
Stp=[StildePrime(i,ron_AA,rcut_AA) for i in t]
plt.plot(t, St, 'r--',t,Stp,'rs',label='AA')
#AB
t = np.arange(ron_AB+0.0001, rcut_AB-0.0001, (rcut_AB-ron_AB)/500)
St=[Stilde(i,ron_AB,rcut_AB) for i in t]
Stp=[StildePrime(i,ron_AB,rcut_AB) for i in t]
plt.plot(t, St, 'b--',t,Stp,'bs',label='AB')
#BB
t = np.arange(ron_BB+0.0001, rcut_BB-0.0001, (rcut_BB-ron_BB)/500)
St=[Stilde(i,ron_BB,rcut_BB) for i in t]
Stp=[StildePrime(i,ron_BB,rcut_BB) for i in t]
plt.plot(t, St, 'g--',t,Stp,'gs',label='BB')
plt.grid()
plt.legend()
plt.show()

#
#Now the full S(r)
#
t = np.arange(0.0001, 1.2*rcut_AA, (2*rcut_AA-0.0001)/500)
St=[S(i,ron_AA,rcut_AA) for i in t]
# Stp=[StildePrime(i,ron_AA,rcut_AA) for i in t]
plt.plot(t, St, 'rs',label='AA')
plt.grid()
plt.show()

#
#The unsmoothened LJ potential
#
plt.figure(2)
plt.xlabel('$r$')
plt.ylabel('$V(r) ~~or~~ V\'(r)$')
#AA
t = np.arange(0.0001, 1.2*rcut_AA, (2*rcut_AA-0.0001)/500)
V=[Vpair(i, eps_AA, sig_AA, rcut_AA) for i in t]
Vp=[VpairPrime(i, eps_AA, sig_AA, rcut_AA) for i in t]
plt.plot(t, V, 'rs',label='$V_{AA}(r)$')
plt.plot(t, Vp, 'r-',label='$V_{AA}\'(r)$')
#AB
t = np.arange(0.0001, 1.2*rcut_AB, (2*rcut_AB-0.0001)/500)
V=[Vpair(i, eps_AB, sig_AB, rcut_AB) for i in t]
Vp=[VpairPrime(i, eps_AB, sig_AB, rcut_AB) for i in t]
plt.plot(t, V, 'bs',label='$V_{AB}(r)$')
plt.plot(t, Vp, 'b-',label='$V_{AB}\'(r)$')
#BB
t = np.arange(0.0001, 1.2*rcut_BB, (2*rcut_BB-0.0001)/500)
V=[Vpair(i, eps_BB, sig_BB, rcut_BB) for i in t]
Vp=[VpairPrime(i, eps_BB, sig_BB, rcut_BB) for i in t]
plt.plot(t, V, 'gs',label='$V_{BB}(r)$')
plt.plot(t, Vp, 'g-',label='$V_{BB}\'(r)$')
axes = plt.gca()
axes.set_ylim([-1.5,5])
plt.legend()
plt.grid()
plt.show()

#
#The unsmoothened LJ potential
#
plt.figure(3)
plt.xlabel('$r$')
plt.ylabel('$V(r) ~~or~~ V\'(r)$')
#AA
t = np.arange(0.0001, 1.2*rcut_AA, (2*rcut_AA-0.0001)/500)
V=[Vij(i, eps_AA, sig_AA, ron_AA, rcut_AA) for i in t]
Vp=[Vprime(i, eps_AA, sig_AA, ron_AA, rcut_AA) for i in t]
plt.plot(t, V, 'rs',label='$V_{AA}(r)$')
plt.plot(t, Vp, 'r-',label='$V_{AA}\'(r)$')
#AB
t = np.arange(0.0001, 1.2*rcut_AB, (2*rcut_AB-0.0001)/500)
V=[Vij(i, eps_AB, sig_AB, ron_AB, rcut_AB) for i in t]
Vp=[Vprime(i, eps_AB, sig_AB, ron_AB, rcut_AB) for i in t]
plt.plot(t, V, 'bs',label='$V_{AB}(r)$')
plt.plot(t, Vp, 'b-',label='$V_{AB}\'(r)$')
#BB
t = np.arange(0.0001, 1.2*rcut_BB, (2*rcut_BB-0.0001)/500)
V=[Vij(i, eps_BB, sig_BB, ron_BB, rcut_BB) for i in t]
Vp=[Vprime(i, eps_BB, sig_BB, ron_BB, rcut_BB) for i in t]
plt.plot(t, V, 'gs',label='$V_{BB}(r)$')
plt.plot(t, Vp, 'g-',label='$V_{BB}\'(r)$')
axes = plt.gca()
axes.set_ylim([-1.5,10])
plt.legend()
plt.grid()
plt.show()

#
#Compare smoothend with unsmoothened
#
plt.figure(2)
plt.xlabel('$r$')
plt.ylabel('$V(r) ~~or~~ V\'(r)$')
#AA
t = np.arange(0.0001, 1.2*rcut_AA, (2*rcut_AA-0.0001)/500)
V=[Vpair(i, eps_AA, sig_AA, rcut_AA) for i in t]
Vp=[VpairPrime(i, eps_AA, sig_AA, rcut_AA) for i in t]
plt.plot(t, V, 'rs',label='$V_{AA}(r)$ NOSHIFT')
plt.plot(t, Vp, 'r-',label='$V_{AA}\'(r)$ NOSHIFT')
t = np.arange(0.0001, 1.2*rcut_AA, (2*rcut_AA-0.0001)/500)
V=[Vij(i, eps_AA, sig_AA, ron_AA, rcut_AA) for i in t]
Vp=[Vprime(i, eps_AA, sig_AA, ron_AA, rcut_AA) for i in t]
plt.plot(t, V, 'bs',label='$V_{AA}(r)$ XPLOR')
plt.plot(t, Vp, 'b-',label='$V_{AA}\'(r)$ XPLOR')
axes = plt.gca()
axes.set_ylim([-1.5,7])
# axes.set_xlim([0.7,1.5])
plt.legend(loc=2)
plt.grid()
plt.show()


raise SystemExit

#Per calcolare la forza dall'integrator, posso sia prenderla dal LJpair, che dalla feature della particella
#Per esempio i dati delle forze di LJ sulla particella 0 si ottengono via:
myLJpair.forces[0] #Non so come si spacchetta per ottenere solo la forza
p = system.particles[0]
p.net_force

#Nello snapshot non vedo la possibilita` di prendere le forze, ma se
#le particelle hanno massa unitaria (come nel mio caso), allora posso usare
acc=snapshot.particles.acceleration

