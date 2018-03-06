#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file at temperature T_ini,
# and: 
# - runs a few steps, to pretend a crappy thermalization.
# - runs until it is clear that at least two ISs have been sampled.
# - takes 2 inherent structures (IS).
# - finds two very close thermal configurations, each leading to a different IS.
# - minimizes each one of them, printing the distance between the two as a 
#   function of the energy of the linear interpolation between the two.
# - some other observables are plottable but not plotted.
#
# To display help:
# python T15-Ridge.py --user="-h"
#
# To launch a simulation:
# python T15-Ridge.py --user="filename -N Natoms --dt=dt <and more arguments>"
# 
# Example:
# python T15-Ridge.py --user="./sample-states/N65.gsd -N65 -T10.0"
#
# List of arguments:
#                filename
# -N             Natoms
# -s             seed (<0: /dev/urandom)
# -t             Number of MD steps done trying to find two different ISs, before giving up.
# -T             temperature
# --tauT         tau of the thermostat
# -d             MD integration step dt
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
from os import remove #Used to remove the backup file at the end of the run
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
from lib import module_measurements as med
from lib import module_timelists as tl
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #
import matplotlib.pyplot as plt

################################################################
#
# SET UP THE SIMULATION CONTEXT
# 
################################################################
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', #positional argument
					nargs=1,
					help='name of the .gsd configuration we want to read'
)
parser.add_argument('-N','--Natoms', #mandatory
					nargs=1,
					type=int,
					required=True,
					help='number of atoms in the system'
)
parser.add_argument('-l','--label', #optional argument
					nargs=1,
					required=False,
					default=['thermalized'],
					help='basename for the output files'
)
parser.add_argument('-s','--seed', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[-1],
					help='seed for random numbers. If negative we get it from /dev/urandom'
)
parser.add_argument('-t','--nNVTsteps', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[10000],
					help='Number of MD steps done trying to find two different ISs, before giving up.'
)
parser.add_argument('-T','--temperature', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[0.6],
					help='target Temperature'
)
parser.add_argument('--tauT', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[1.0],
					help='tau of the thermostat'
)
parser.add_argument('--dt', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[0.0025],
					help='dt for MD integration'
)


args = parser.parse_args(more_arguments)
filename=args.filename[0]
Natoms=args.Natoms[0]
label=args.label[0]
seed=args.seed[0]
nNVTsteps=args.nNVTsteps[0]
temperature=args.temperature[0]
tauT=args.tauT[0]
dt=args.dt[0]
dtFIRE=dt/10
THRES=1e-8
del parser

if seed>0:
   np.random.seed(seed)
print("Input configuration: ",filename)
print("Natoms = ",Natoms)
print("seed = ",seed)
print("nNVTsteps = ",nNVTsteps)
print("T = ",temperature)
print("tauT = ",tauT)
print("dt = ",dt)
print("dtFIRE = ",dtFIRE)
print("thermostat = NVT")
print("label = ",label)
assert(nNVTsteps>0)
assert(temperature>=0)
assert(tauT>0)
assert(dt>0 and dt<0.1)

################################################################
#
# CLASS AND FUNCTION DEFINITION
#
################################################################
def SetupAnalyzer(logname="log.txt", period=2000, quantities=['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'],seed=seed):
	analyzer = hoomd.analyze.log(filename=logname, \
								 quantities=quantities, period=period, \
								 header_prefix = '#seed:'+str(seed)+"\n#", \
								 overwrite=False,
								 phase=0)
	return analyzer

def PotEn(mode,analyzer,dt=0.0025):
	mode.dt=1e-18
	hoomd.run(2)
	U=analyzer.query('potential_energy')
	mode.dt=dt
	return U

def Minimize(snap,verbose=True):
	with simIS:
		systemIS.restore_snapshot(snap)
		modeFire.reset()
		while not(modeFire.has_converged()):
			hoomd.run(100)
		eIS=PotEn(modeFire,analyzerFire,dt=dtFIRE)
		if verbose: print("EIS(T) = ",eIS)
		snapIS=systemIS.take_snapshot(dtype='double')
	return snapIS,eIS

def EvolveOneStep(verbose=True):
	with simT:
		hoomd.run(1)
		eT=PotEn(modeStandard,analyzerStandard,dt=dt)
		if verbose:	print("t: ",t, "E(T) = ",eT, end='\t')
		snapT=systemT.take_snapshot(dtype='double')
	return snapT,eT

def LinearConfInterpolation(snap1, snap2, box_size):
	"""
	returns a linear interpolation between snap1 and snap2
	new_positions = (pos1 + pos2)/2
	"""
	snap12=systemIS.take_snapshot()
	snap12.particles.position[:] = med.PeriodicIntermPoints(
		np.array(snap1.particles.position,dtype=np.float64),
		np.array(snap2.particles.position,dtype=np.float64),
		box_size)
	return snap12

def ConfBisect(snap1, snap2, eis1, eis2, L, dmax=0.002):
	"""
	This function takes two thermal snapshots that end in two separate IS,
	and makes a linear interpolation between the positions of the two.
	Through bisection, the point that separates the two basins is found.
	"""
	assert(np.abs(eis1-eis2)>THRES)
	dstart=0.1*dmax
	Natoms=snap1.particles.N
	dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms #the box is cubic
	snap12=LinearConfInterpolation(snap1, snap2, L)
	snapis12,eis12=Minimize(snap12,verbose=False)

	count=0
	maxcount=100
	while dist12>dstart:
		if np.abs(eis1-eis12) <= THRES: #If snap12 belongs to snap1, snap1=snap12
			snap1.particles.position[:]=snap12.particles.position
			eis1=eis12
			pos1=np.array(snap1.particles.position[:], dtype=np.float64)
		elif np.abs(eis2-eis12) <= THRES: #If snap12 belongs to snap2, snap2=snap12
			snap2.particles.position[:]=snap12.particles.position
			eis2=eis12

		else: #If snap12 does not belong to either, we throw a warning and change snap2
			print("ConfBisect: found an intermediate IS while searching the TS")
			snap2.particles.position[:]=snap12.particles.position
			eis2=eis12

		dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms
		count+=1
		assert(np.abs(eis1-eis2)>THRES)
		snap12=LinearConfInterpolation(snap1, snap2, L)
		snapis12,eis12=Minimize(snap12,verbose=False)

		if count>maxcount:
			raise SystemExit("ConfBisect ERROR: the interpolation bisection is not converging.")
	return snap1,snap2,snap12,eis1,eis2,eis12,dist12


################################################################
#
# INITIALIZE
#
################################################################
backupname=label+"_backup.gsd"
simT  = hoomd.context.SimulationContext();
simIS = hoomd.context.SimulationContext();

with simT:
	systemT = hoomd.init.read_gsd(filename=filename, restart=None)
	assert(Natoms==len(systemT.particles))
	myLjPair=pot.LJ(md.nlist.cell(),type="KAshort")
	analyzerStandard=SetupAnalyzer(logname=None, period=None)
	modeStandard=md.integrate.mode_standard(dt=dt)
	md.update.zero_momentum(phase=2, period=10)
	integrator = md.integrate.nvt(group=hoomd.group.all(), kT=temperature, tau=tauT)
	hoomd.run(5/dt)
	snapT=systemT.take_snapshot(dtype='double')
	boxParamsT=snapT.box

with simIS:
	systemIS = hoomd.init.read_gsd(filename=filename, restart=None)
	assert(Natoms==len(systemIS.particles))
	myLjPair=pot.LJ(md.nlist.cell(),type="KAshort")
	analyzerFire=SetupAnalyzer(logname=None, period='None')
	modeFire=md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=0.99, ftol=1e-5, Etol=1e-12, wtol=1e-5)
	md.update.zero_momentum(phase=2, period=10)
	integrator = md.integrate.nve(group=hoomd.group.all())

L=np.float64(boxParamsT.Lx)
assert(np.abs(boxParamsT.Lx-boxParamsT.Lz)<THRES)

################################################################
# 
# INTEGRATION
# 
################################################################

################################################
# Run dynamics until two separate IS are found #
################################################

print("-- Run dynamics until two distinct ISs are found --")

eISold=0
for t in range(nNVTsteps):
	snapT,eT=EvolveOneStep()
	snapIS,eIS=Minimize(snapT)

	if np.abs(eISold-eIS)>THRES and t>0:
		break
	(eISold,eTold,snapISold,snapTold)=(eIS,eT,snapIS,snapT)
if t==nNVTsteps-1: sys.exit("Unable to find two distinct inherent structures")

Eis1=eISold
Eis2=eIS
ET1=eTold
ET2=eT
snapis1=snapISold
snapis2=snapIS
snapT1=snapTold
snapT2=snapT
snapT1.particles.velocity[:] = np.zeros((Natoms, 3))
snapT2.particles.velocity[:] = np.zeros((Natoms, 3))


print("ET1 = ",ET1)
print("ET2 = ",ET2)
print("Eis1 = ",Eis1)
print("Eis2 = ",Eis2)

###################
#                 #
# Calculate Ridge #
#                 #
###################
print("-- Find two thermal configurations that are very close to eachother, but end up in two different ISs --")
dmax=0.0004 #0.004 is about the typical distance between confs at subsequent time steps with dt=0.0025
snap1,snap2,snap12,eis1,eis2,eis12,dist12=ConfBisect(snapT1, snapT2, Eis1, Eis2, L, dmax=dmax)
dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms
print("eis1 = ", eis1,"\teis2 = ",eis2)



print("-- Prepare two parallel minimization contexts --")
snapIS1=snap1
snapIS2=snap2
snapIS1.particles.velocity[:]=np.zeros((Natoms, 3))
snapIS2.particles.velocity[:]=np.zeros((Natoms, 3))
simIS1 = hoomd.context.SimulationContext();
simIS2 = hoomd.context.SimulationContext();
with simIS1:
	systemIS1 = hoomd.init.read_snapshot(snapIS1)
	assert(Natoms==len(systemIS1.particles))
	myLjPair=pot.LJ(md.nlist.cell(),type="KAshort")
	analyzerFire1=SetupAnalyzer(logname=None, period='None')
	modeFire1=md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=0.99, ftol=1e-5, Etol=1e-12, wtol=1e-5)
	md.update.zero_momentum(phase=2, period=10)
	integrator1 = md.integrate.nve(group=hoomd.group.all())
with simIS2:
	systemIS2 = hoomd.init.read_snapshot(snapIS2)
	assert(Natoms==len(systemIS2.particles))
	myLjPair=pot.LJ(md.nlist.cell(),type="KAshort")
	analyzerFire2=SetupAnalyzer(logname=None, period='None')
	modeFire2=md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=0.99, ftol=1e-5, Etol=1e-12, wtol=1e-5)
	modeFire2.reset()
	md.update.zero_momentum(phase=2, period=10)
	integrator2 = md.integrate.nve(group=hoomd.group.all())




print("-- Minimize the two configurations in parallel --")
nsteps=10
niter=10000

distances12=[dist12]
distances1=[med.PeriodicDistance(snap1.particles.position, snapIS1.particles.position, L).sum()/Natoms]
distances2=[med.PeriodicDistance(snap2.particles.position, snapIS2.particles.position, L).sum()/Natoms]
energies1=[ET1]
energies2=[ET2]
energies12=[myLjPair.CalculateEnergySlower(LinearConfInterpolation(snapIS1, snapIS2, L))]
steps=[0]
modeFire1.reset()
modeFire2.reset()
eRidgeFound=False
for iter in range(niter):
	with simIS1:
		systemIS1.restore_snapshot(snapIS1)
		hoomd.run(nsteps)
		eIS1=PotEn(modeFire1,analyzerFire1,dt=dtFIRE)
		snapIS1=systemIS1.take_snapshot(dtype='double')
	print((iter+1)*nsteps,"\teIS1 = ",eIS1,end='\t')
	with simIS2:
		systemIS2.restore_snapshot(snapIS2)
		hoomd.run(nsteps)
		eIS2=PotEn(modeFire2,analyzerFire2,dt=dtFIRE)
		snapIS2=systemIS2.take_snapshot(dtype='double')
	print("eIS2 = ",eIS2,end='\t')
	dist12=med.PeriodicDistance(snapIS1.particles.position, snapIS2.particles.position, L).sum()/Natoms
	distances12.append(dist12)
	distances1.append(med.PeriodicDistance(snap1.particles.position, snapIS1.particles.position, L).sum()/Natoms)
	distances2.append(med.PeriodicDistance(snap2.particles.position, snapIS2.particles.position, L).sum()/Natoms)
	energies1.append(eIS1)
	energies2.append(eIS2)
	energy12=myLjPair.CalculateEnergySlower(LinearConfInterpolation(snapIS1, snapIS2, L))
	energies12.append(energy12)
	steps.append((iter+1)*nsteps)
	print("dist = ",dist12)
	if (not eRidgeFound) and dist12>0.001:
		Eridge=energy12
		eRidgeFound=True #Se l'obiettivo fosse solo trovare Eridge, invece di fare un tutorial, metteremmo un break qui.
	if modeFire1.has_converged() and modeFire2.has_converged(): break

if iter==niter-1: sys.exit("Unable to find Eridge (minimizations did not converge). Try with niter>>",niter)

print("--Finished!--")


print("-- Plots --")
print("We found Er = ",Eridge)
print("Let us see if the graphs confirm it")

#Invert the signs of the energies, to be able to use log scales
energies1=-np.array(energies1)
energies2=-np.array(energies2)
energies12=-np.array(energies12)
Er=-Er



plt.semilogy(energies12,distances12,label='distance (1 from 2)')
axes = plt.gca()
axes.set_xscale("linear", nonposx='clip')
axes.set_yscale("linear", nonposy='clip')
plt.title('$T = $'+str(temperature))
plt.xlabel('$E_{12}$')
plt.ylabel('$d_{12}$')
plt.grid(True)
plt.legend()
axes.set_xlim(Eridge-1,Eridge+1)
plt.savefig("dist-E_T"+str(temperature)+".png")
plt.show()
