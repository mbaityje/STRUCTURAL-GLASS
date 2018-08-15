#!/usr/bin/env python
################################################################
#
#
# DESCRIPTION
# Given two minima, finds the path between them through the Nudged Elastic Band method.
# This is a very basin implementation.
# Backtracking is implemented, in the form suggested in J. Chem. Theory Comput., 2017, 13 (7), pp 3250–3259
# 
# The algorithm is slow: launch on small systems.
# 
# To display help:
# python T27-NudgedElasticBandBacktrack.py --user="-h"
#
# To launch a simulation:
# python T27-NudgedElasticBandBacktrack.py --user="conf1 conf2 --f1=frame1 --f2=frame2 -n=npoints --maxiter=maxiterNEB --alpha=initalpha -k=springConstant"
#
# For example:
# python T27-NudgedElasticBandBacktrack.py --user="./sample-states/rotenbergKA_T2.0_N1080.gsd ./sample-states/rotenbergKA_T2.0_N1080.gsd -n20"
# python T27-NudgedElasticBandBacktrack.py --user="../OUTPUT/T0.6/N65/S3/restartChunk139.gsd ../OUTPUT/T0.6/N65/S3/restartChunk141.gsd -n10 --maxiter=200 -k1 --alpha=1e-5"
################################################################


#
# LIBRARIES, OPTIONS AND HARD-CODED PARAMETERS
#
from __future__ import print_function 
import sys, argparse, copy, math
from matplotlib import pyplot as plt
import numpy as np 
import hoomd, gsd.pygsd, gsd.hoomd
from lib import module_potentials as pot, module_measurements as med
np.set_printoptions(precision=15)
fireParams={'alpha':0.99, 'ftol':1e-5, 'Etol':1e-10, 'wtol':1e-5, 'minsteps':100}

#
# FUNCTIONS
#
def Init():
	'''
	HOOMD INITIALIZATION AND COMMAND-LINE ARGUMENTS	
	'''
	hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
	hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
	more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd
	parser = argparse.ArgumentParser(add_help=True)
	parser.add_argument('filename1', help='name of the first configuration file .gsd')
	parser.add_argument('filename2', help='name of the second configuration file .gsd')
	parser.add_argument('--f1', type=int, required=False, default=0, help='read the f1-th frame of filename1')
	parser.add_argument('--f2', type=int, required=False, default=0, help='read the f2-th frame of filename2')
	parser.add_argument('-n','--npoints', type=int, required=False, default=10, help='number of points in the interpolation')
	parser.add_argument('--maxiter', type=int, required=False, default=10, help='maximum NEB iterations')
	parser.add_argument('--alpha', type=float, required=False, default=-1.0, help='(initial) value of the learning rate')
	parser.add_argument('-k', type=float, required=False, default=1, help='spring constant for the NEB')
	args = parser.parse_args(more_arguments)
	return args

def ReadConfs(args):
	system=hoomd.init.read_gsd(args.filename1, frame=args.f1, time_step=0, restart=None)
	snap=system.take_snapshot()
	Natoms = snap.particles.N
	L=snap.box.Lx

	#Read Second Configuration
	with open(args.filename2, 'rb') as flow:
		HoomdFlow = gsd.pygsd.GSDFile(flow)
		hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
		s2=hoomdTraj.read_frame(args.f2) 

	# Check that the two configurations represent the same system
	# Same number of particles
	if s2.particles.N != Natoms: 
		raise ValueError('The two configurations have different numbers of particles (%d and %d)'%(Natoms,s2.particles.N))
	#Same particle ids
	if all(snap.particles.typeid==s2.particles.typeid) == False: 
		raise ValueError('The two configurations do not have the same particle indices')
	#Same space dimensionality
	if snap.box.dimensions != s2.configuration.dimensions: 
		raise ValueError('The two configurations live in different spatial dimensions (%d and %d)'%(snap.box.dimensions, s2.configuration.dimensions))
	#Same box size
	if not (math.isclose(snap.box.Lz,L) and math.isclose(snap.box.Ly,L) ): 
		raise ValueError('The box must be a square because of the implementation of periodic boundary conditions.') 
	if L==s2.configuration.box[0]==False: 
		raise ValueError('Either the two boxes are different (Lx1=%g and Lx2=%g), or you should make a decent floating point comparison'%(snap.box.Lx, s2.configuration.box[0]))

	return system, Natoms, np.array(snap.particles.position[:], dtype=np.float64), np.array(s2.particles.position[:], dtype=np.float64), L


def InitSystem(Natoms):
	NeighborsList = hoomd.md.nlist.cell()
	mypot="KAshort" if Natoms<500 else "KA"
	potential=pot.LJ(NeighborsList,type=mypot)
	analyzer = hoomd.analyze.log(filename=None, quantities=['temperature','potential_energy', 'kinetic_energy', 'momentum'], period=None, header_prefix = '#', overwrite=True, phase=0)
	return potential,analyzer

def Minimize(snap):
	system.restore_snapshot(snap)
	fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=fireParams['alpha'], ftol=fireParams['ftol'], Etol=fireParams['Etol'], wtol=fireParams['wtol'], min_steps=fireParams['minsteps'])
	integrator=hoomd.md.integrate.nve(group=hoomd.group.all())
	fire.cpp_integrator.reset()
	while not(fire.has_converged()):
		hoomd.run(100)
	integrator.disable()
	return analyzer.query('potential_energy')

def MinimizeConfs(system, pos1, pos2):

	def MinimizeOneConf(system, pos):
		snap=system.take_snapshot()
		snap.particles.velocity[:]=np.zeros((Natoms,3))
		snap.particles.position[:]=pos[:]
		system.restore_snapshot(snap)
		eis=Minimize(snap)
		snap=system.take_snapshot()
		print("Eis=",eis)
		return eis, np.array([system.particles[i].position for i in range(len(system.particles))], dtype=np.float64)

	E1,pos1=MinimizeOneConf(system, pos1)
	E2,pos2=MinimizeOneConf(system, pos2)

	return E1, E2, pos1, pos2

def InterpolatePivots(pos1, pos2, npoints, L):
	'''
	Interpolate npoints between pos1 and pos2. L is the box size.
	'''
	delta=1./npoints
	return np.array([med.PeriodicInterpPoints(pos2, pos1, L, i*delta) for i in range(npoints+1)], dtype=np.float64)

def PathEnergies(mypivots):
	energies=[]
	snap=system.take_snapshot()
	for i in range(len(mypivots)):
		snap.particles.position[:]=mypivots[i][:]
		energy=potential.CalculateEnergySlower(snap)
		# energy=EnergyFromAnalyzer()
		energies.append(energy)
	return np.array(energies)

def EnergyFromAnalyzer():
	hoomd.run(2)
	return analyzer.query('potential_energy')

def EnergyAccFromAnalyzer():
	hoomd.run(2)
	return analyzer.query('potential_energy'), np.array([system.particles[i].acceleration for i in range(len(system.particles))], dtype=np.float64)

def MeasureDistances(pivots):
	'''
	Misura le distanze tra i pivot in ascissa curvilinea
	'''
	distances=[0]
	for i in range(1,len(pivots)):
		distances.append(distances[i-1]+med.PeriodicDistance(pivots[i], pivots[i-1], L) )
	return np.array(distances, dtype=np.float64)


def CalcTangent(i, pivots, Eprev, Ethis, Enext, L):
	dispPrev=med.PeriodicDisplacement(pivots[i],pivots[i-1],L)
	dispNext=med.PeriodicDisplacement(pivots[i+1],pivots[i],L)
	dispPrevNorm=np.linalg.norm(dispPrev)
	dispNextNorm=np.linalg.norm(dispNext)

	if Enext>Ethis != (Eprev>Ethis): #is not an extremum
		tangent = dispNext/dispNextNorm if Enext>Eprev else dispPrev/dispPrevNorm
	else: #if it is an extremum, apply smoothening
		delta=np.abs([Enext-Ethis, Eprev-Ethis])
		deltamax=np.max(delta)
		deltamin=np.min(delta)
		tangent = dispNext*deltamax+dispPrevNorm*deltamin if Enext>Eprev else dispNext*deltamin+dispPrevNorm*deltamax
		tangent/= np.linalg.norm(tangent)
	return tangent, dispPrevNorm, dispNextNorm



def BackTrack(fRMS, fRMSold, myalpha, mynBack, fmax, L, eps=1e-3, gamma=0.9):
	assert(gamma<1)
	n0=5
	chk=(fRMS-fRMSold)/np.abs(fRMS+fRMSold)
	if chk>eps:
		myalpha*=gamma
		mynBack=n0
	else:
		mynBack-=1
		if mynBack<0:
			myalpha=min(myalpha/gamma, L/fmax*0.5) #No alpha should displace a particle more than half the box size
			myalpha/=gamma
			mynBack=n0

	return myalpha, mynBack

# ###################### #
# ###################### #
# HERE STARTS THE SCRIPT #
# ###################### #
# ###################### #

# Initialization
args = Init()
system,Natoms,pos1,pos2,L = ReadConfs(args)
potential,analyzer = InitSystem(Natoms)
Eini, Efin, pos1, pos2 = MinimizeConfs(system, pos1, pos2)


#
# CALCULATE INITIAL INTERPOLATION (GUESS OF THE PATH)
#
mode=hoomd.md.integrate.mode_standard(dt=1e-20)
integrator=hoomd.md.integrate.nve(group=hoomd.group.all())
pivots=InterpolatePivots(pos1, pos2, args.npoints+2, L)
energiesOld=PathEnergies(pivots)
distances=MeasureDistances(pivots)
ax = plt.axes()
ax.set_color_cycle([plt.cm.cool(t) for t in np.linspace(0, 1, args.maxiter)])
plt.plot(distances, energiesOld, marker='.',linestyle='-', label='$t$ = 0')
plt.plot(distances, energiesOld, label='$t$ = 0')


#
# NUDGED ELASTIC BAND
#
#Set dummy integrator for measurements
snap=system.take_snapshot()
snap.particles.velocity[:]=np.zeros((Natoms,3))

alpha=np.full(len(pivots), args.alpha, dtype=np.float32)
nBack=np.full(len(pivots), 2, dtype=int)
fRMSold=np.ndarray(len(pivots))
for t in range(1,args.maxiter):
	#Redistribute confs?

	print("t: ",t, end='\t')

	Eprev=Eini
	snap.particles.position[:]=pivots[1][:]
	system.restore_snapshot(snap)
	Ethis,Athis=EnergyAccFromAnalyzer()

	error=0
	for i in range(1, len(pivots)-1):
		#CalcNextEnergy
		snap.particles.position[:]=pivots[i+1][:]
		system.restore_snapshot(snap)
		Enext,Anext=EnergyAccFromAnalyzer()
		#CalcTangent
		tangent,dispPrevNorm,dispNextNorm=CalcTangent(i, pivots, Eprev, Ethis, Enext, L)
		#CalcNEBForce
		confForcePerp=Athis-np.sum(Athis*tangent)*tangent #F=A because m=1
		springForce=args.k*(dispNextNorm-dispPrevNorm)*tangent
		#Update
		force=confForcePerp+springForce
		fRMS=np.linalg.norm(force)
		fmax=np.abs(force).max()
		if t==1:
			if args.alpha<0: #If the user didnt set the initial alpha, we choose it now
				alpha[i]=L/fmax*0.001 #In the first iteration the maximum displacement is 0.1% of the box
		else:
			alpha[i], nBack[i] = BackTrack(fRMS, fRMSold[i], alpha[i], fmax, L, nBack[i])
		pivots[i] = med.PeriodicSum(pivots[i], alpha[i]*force, L)

		#Prepare for next iteration
		Eprev=Ethis
		Ethis=Enext
		Aprev=Athis
		Athis=Anext
		fRMSold[i]=fRMS

	distances=MeasureDistances(pivots)
	energies=PathEnergies(pivots)
	error=np.abs(energies-energiesOld).sum()/args.npoints
	print('error: ',error)
	if error<1e-3:
		break
	energiesOld=energies
	plt.plot(distances, energies,'.', label='$t$ = '+str(t))

plt.plot(distances, energies, label='$t$ = '+str(t))
plt.xlabel('distance')
plt.ylabel('energy')
# plt.legend()
plt.savefig('./test-output/neb_all.png')
plt.show()

plt.plot(distances, energies, label='$t$ = '+str(t))
plt.xlabel('distance')
plt.ylabel('energy')
# plt.legend()
plt.savefig('./test-output/neb_final.png')
plt.show()

