#!/usr/bin/env python3
################################################################
#
# DESCRIPTION
# 
# Reads a thermal trajectory starting from frame $iframe.
# Minimizes the following $nframes frames.
# Calculates Eis(jframe), msdIS(jframe-1,jframe), nIS(jframe-1,jframe) [# of particles that moved in time dt].
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import argparse
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import module_measurements as med
import os.path
import gsd.pygsd
import gsd.hoomd

################################################################
#
# READ ARGUMENTS
# 
################################################################
#Start hoomd
print("Initialize hoomd context\n")
simT=hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user()

#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', nargs=1, help='name of the initial state .gsd')
parser.add_argument('--nframes', nargs=1, type=int, required=True, help='Total number of frames')
parser.add_argument('--iframe', nargs=1, type=int, required=False, default=[0], help='Frame to read from the gsd file (default is 0)')
parser.add_argument('--saveISgsd', nargs=1, type=bool, required=False, default=[False], help='If True, saves the IS trajectory')
parser.add_argument('--saveThermalgsd', nargs=1, type=bool, required=False, default=[False], help='If True, saves the thermal trajectory')
parser.add_argument('-l','--label', nargs=1, required=False, default=[''], help='label for distinguishing runs and continuations')
args = parser.parse_args(more_arguments)

filename=args.filename[0]
nframes=args.nframes[0]
iframe=args.iframe[0]
saveISgsd=args.saveISgsd[0]
saveThermalgsd=args.saveThermalgsd[0]
label=str(args.label[0])
THRES=1e-6

#Integrator parameters
dtFIRE=0.0025
alphaFIRE=0.99
ftolFIRE=1e-5
EtolFIRE=1e-10
wtolFIRE=1e-5
minstepsFIRE=100

print("filename: ",filename)
print("initial frame = ",iframe)
print("nframes = ",nframes)
if nframes <0: raise ValueError('Received nframes='+str(nframes)+'nframes must be >=0')



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
	print("Natoms = ",Natoms)
	boxParams=s0.configuration.box
	L=boxParams[0]
	if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
		print('box dimensions are : ', boxParams[0])
		print('and are not supported (only isotropic systems supported)')
		raise SystemExit

	#Make sure that we don't ask for more frames than those present in the trajectory    
	totframes = len(hoomdTraj)
	if(iframe+nframes>totframes):
		nframes=totframes-iframe
		finalframe=totframes
		print("Shortening nframes. Now nframes=",nframes)
	else:
		finalframe=iframe+nframes
	posizioni=[hoomdTraj[i].particles.position[:] for i in range(iframe,finalframe)]
	HoomdFlow.close()


################################################################
# 
# Initialize
#
################################################################
system = hoomd.init.read_gsd(filename=filename)
snapT_old=system.take_snapshot()
snapT_old.particles.position[:]=posizioni[0]
snapT_old.particles.velocity[:] = np.zeros((Natoms,3))
snapT=system.take_snapshot()
snapT.particles.velocity[:] = np.zeros((Natoms,3))

################################################################
# 
# Set potential
#
################################################################
NeighborsListLJ = md.nlist.cell()
potential=pot.LJ(NeighborsListLJ,type="KAshort")

################################################################
# 
# Set analyzer
#
################################################################
analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
analyzer = hoomd.analyze.log(filename=None, quantities=analyzer_quantities, period=None)


################################################################
# 
# Functions
#
################################################################

def Minimize(snap):
	system.restore_snapshot(snap)
	fire.cpp_integrator.reset()
	if not integrator_fire.enabled: integrator_fire.enable()
	while not(fire.has_converged()):
		hoomd.run(100)
	eIS=analyzer.query('potential_energy')
	snapnew=system.take_snapshot()
	integrator_fire.disable()
	return eIS,snapnew

def LinearConfInterpolation(snap1, snap2, box_size):
	"""
	returns a linear interpolation between snap1 and snap2
	new_positions = (pos1 + pos2)/2
	"""
	snap12=system.take_snapshot()
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
	#0.002 is about half the typical distance between confs at subsequent time steps w/ dt=0.0025
	"""
	print('···ConfBisect···')
	print('eis1 = ',eis1,'; eis2 = ',eis2)

	if np.abs(eis1-eis2)<THRES: raise SystemExit('ConfBisect ERROR: the two starting ISs the same one [abs(eis1-eis2)<THRES]')
	dstart=0.5*dmax
	Natoms=snap1.particles.N
	dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms #the box is cubic
	snap12=LinearConfInterpolation(snap1, snap2, L)

	eis12,snapis12=Minimize(snap12)

	count=0
	maxcount=100
	while dist12>dstart:
		print('eis1 = ',eis1,'; eis2 = ',eis2,'; eis12 = ',eis12)
		if np.abs(eis1-eis12) <= THRES: #If snap12 belongs to snap1, snap1=snap12
			snap1.particles.position[:]=snap12.particles.position
			eis1=eis12
		elif np.abs(eis2-eis12) <= THRES: #If snap12 belongs to snap2, snap2=snap12
			snap2.particles.position[:]=snap12.particles.position
			eis2=eis12
		else: #If snap12 does not belong to either, we print a warning and change snap2
			print("ConfBisect: found an intermediate IS while searching the TS (Eis=%.14g)"%eis12)
			snap2.particles.position[:]=snap12.particles.position
			eis2=eis12
		dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms
		count+=1
		if  np.abs(eis1-eis2)<THRES: raise SystemExit('ConfBisect ERROR: the two ISs became the same one [abs(eis1-eis2)<THRES]')
		snap12=LinearConfInterpolation(snap1, snap2, L)
		eis12,snapis12=Minimize(snap12)

		if count>maxcount:
			raise RecursionError("ConfBisect ERROR: the interpolation bisection is not converging.")
	return snap1,snap2,snap12,eis1,eis2,eis12,dist12

def SetupAnalyzer(logname=None, period=None, quantities=['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum']):
	analyzer = hoomd.analyze.log(filename=logname, \
								 quantities=quantities, period=period, \
								 header_prefix = '#', \
								 overwrite=False,
								 phase=0)
	return analyzer



def CalculateRidge(snapT1, snapT2, Eis1, Eis2, L):
	''' Calculates energy at the ridge between two snapshots '''
	print("--- Calculate Ridge ---\n- Eis1:",Eis1,' Eis2:',Eis2)

	dmax=0.0001 #0.004 is about the typical distance between confs at subsequent time steps with dt=0.0025
	nsteps=1
	niter=10000
	snapT1.particles.velocity[:]=np.zeros((Natoms, 3))
	snapT2.particles.velocity[:]=np.zeros((Natoms, 3))

	#snapis are not inherent structures, but the gradually approach them
	snapis1,snapis2,snapis12,Eis1,Eis2,Eis12,dist12=ConfBisect(snapT1, snapT2, Eis1, Eis2, L, dmax=dmax)

	print("- Eis1=%.8f;\tEis2= %.8f (after ConfBisect)"%(Eis1,Eis2))

	context1 = hoomd.context.SimulationContext();
	with context1:
		system1 = hoomd.init.read_snapshot(snapis1)
		modeFire1=md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=alphaFIRE, ftol=ftolFIRE, Etol=EtolFIRE, wtol=wtolFIRE)
		integrator1 = md.integrate.nve(group=hoomd.group.all())
		analyzerFire1=SetupAnalyzer(logname=None, period='None')
		md.update.zero_momentum(phase=2, period=10)
		pair=pot.LJ(md.nlist.cell(),type="KAshort")

	context2 = hoomd.context.SimulationContext();
	with context2:
		system2 = hoomd.init.read_snapshot(snapis2)
		modeFire2=md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=alphaFIRE, ftol=ftolFIRE, Etol=EtolFIRE, wtol=wtolFIRE)
		integrator2 = md.integrate.nve(group=hoomd.group.all())
		analyzerFire2=SetupAnalyzer(logname=None, period='None')
		md.update.zero_momentum(phase=2, period=10)
		pair=pot.LJ(md.nlist.cell(),type="KAshort")

	modeFire1.reset()
	modeFire2.reset()


	for iter in range(niter):
		with context1:
			system1.restore_snapshot(snapis1)
			hoomd.run(nsteps)
			eis1=analyzerFire1.query('potential_energy')
			snapis1=system1.take_snapshot(dtype='double')
		with context2:
			system2.restore_snapshot(snapis2)
			hoomd.run(nsteps)
			eis2=analyzerFire2.query('potential_energy')
			snapis2=system2.take_snapshot(dtype='double')

		dist12=med.PeriodicDistance(snapis1.particles.position, snapis2.particles.position, L).sum()/Natoms

		print("iter: ",iter," dist12=",dist12,"eis1: %.14f"%eis1," eis2: %.14f"%eis2)

		if dist12>dmax:
			'''
			With the linear interpolations the found barrier is occasionally lower than one of the ISs, because of nonlinearities in the path.
			# snapRidge=LinearConfInterpolation(snapis1_old, snapis2_old, L)
			# Eridge=potential.CalculateEnergySlower(snapRidge)
			For this reason, I just take the energy of one of the two confs at the previous step.
			'''
			if eis1>eis2: 
				Eridge=ConsistentRidge(eis1_old,Eis1,Eis2,THRES)
				snapRidge=snapis1_old
			else : 
				Eridge=ConsistentRidge(eis2_old,Eis1,Eis2,THRES)
				snapRidge=snapis2_old
			print("Eridge = ",Eridge)
			break

		snapis1_old=snapis1
		snapis2_old=snapis2
		eis1_old=eis1
		eis2_old=eis2

	if iter==niter-1: sys.exit('CalculateRidge did not converge after '+str(niter)+' steps')
	return Eridge,snapRidge


def ConsistentRidge(Eridge, Eis1, Eis2, thres):
	'''
	If Eridge>Eis1 and Eridge>Eis2 everything is consistent.
	For precision reasons, it might happen that we must use the condition Eridge+thres>Eis.
	If neither this condition is met, something is wrong.
	'''
	if Eridge>Eis1:
		if Eridge>Eis2: 
			return Eridge
		elif Eridge+thres>Eis2:
			return Eridge+thres
		else:
			sys.exit('Inconsistent Eridge<Eis.\nEridge='+str(Eridge)+';\tEis2='+str(Eis2))
	elif Eridge+thres>Eis1:
		return Eridge+thres
	else:
		sys.exit('Inconsistent Eridge<Eis.\nEridge='+str(Eridge)+';\tEis1='+str(Eis1))

def EnerSnapshot(snapshot):
	system.restore_snapshot(snapshot)
	modeT.dt=1e-24
	if not integrator.enabled: integrator.enable()
	hoomd.run(2)
	ener=analyzer.query('potential_energy')
	integrator.disable()
	return ener

# def PotEn(mode,analyzer,dt=0.0025):
# 	mode.dt=1e-24
# 	hoomd.run(2)
# 	U=analyzer.query('potential_energy')
# 	mode.dt=dt
# 	return U



################################################################
# 
# Minimizations
#
################################################################

#Define lists of observables
EISlist=[]
ETlist=[]
EridgeList=[]
tRidgeList=[]
msdlist=[]
nmovedlist=[]

#Integrators
modeT=md.integrate.mode_standard(dt=1e-24)
integrator = md.integrate.nve(group=hoomd.group.all())

#Initial thermal energy
ET=EnerSnapshot(snapT_old)
if saveThermalgsd==True:
	confnameTherm='MinimizeSegmentTherm.gsd'
	hoomd.dump.gsd(confnameTherm, period=None, group=hoomd.group.all(), overwrite=True, truncate=False, phase=-1)

# Initial IS
fire=hoomd.md.integrate.mode_minimize_fire(dt=dtFIRE, alpha_start=alphaFIRE, ftol=ftolFIRE, Etol=EtolFIRE, wtol=wtolFIRE, min_steps=minstepsFIRE)
integrator_fire = md.integrate.nve(group=hoomd.group.all())

Eis_old,snapIS_old=Minimize(snapT_old)
print("ET=",ET," Eis=",Eis_old)


#
# Minimize for further times
#
totmoved=0
for jframe in range(iframe+1,finalframe):
	print("jframe:",jframe,"jframe-iframe:",jframe-iframe)

	#Update configuration
	snapT_old.particles.position[:]=posizioni[jframe-iframe-1]
	snapT.particles.position[:]=posizioni[jframe-iframe]

	#Thermal energy
	ET=EnerSnapshot(snapT)
	ETlist.append(ET)
	if saveThermalgsd==True: hoomd.dump.gsd(confnameTherm, period=None, group=hoomd.group.all(), overwrite=False, truncate=False, phase=-1)

	#Inherent structure
	Eis,snapIS=Minimize(snapT)
	if saveISgsd==True: hoomd.dump.gsd(confnameIS, period=None, group=hoomd.group.all(), overwrite=False, truncate=False, phase=-1)

	if np.abs(Eis - Eis_old)>THRES:
		tRidgeList.append(jframe-0.5) #Ridge is between the two steps
		Eridge,snapRidge=CalculateRidge(snapT_old, snapT, Eis_old, Eis, L)
		print('>>>>>>>>>>>>>>>>')
		print('Eridge = ',Eridge)
		print('>>>>>>>>>>>>>>>>')
		EridgeList.append(Eridge)

	msdIS=med.PeriodicSquareDistance(snapIS_old.particles.position, snapIS.particles.position, L)
	nmoved=med.HowManyMovedPos(snapIS_old.particles.position, snapIS.particles.position,L)
	print("ET=",ET," Eis=",Eis,"  msdIS=",msdIS,"  nmoved=",nmoved)
	totmoved+=nmoved

	EISlist.append(Eis)
	msdlist.append(msdIS)
	nmovedlist.append(nmoved)

	snapIS_old=snapIS
	Eis_old=Eis
print("totmoved=",totmoved)


#Output of IS energies
output=np.column_stack((range(iframe+1,finalframe), ETlist, EISlist, msdlist, nmovedlist))
np.savetxt('MinimizeSegmentWithRidges'+label+'.txt', output,fmt='%d %.14g %.14g %.14g %.14g', header='1)iframe 2)E(T)(iframe) 3)Eis(iframe) 4)msd(iframe-1,iframe) 5)nmoved(iframe-1,iframe)')


#output of ridges
Eante=[EISlist[int(i-0.5)-iframe] for i in tRidgeList]
Epost=[EISlist[int(i+0.5)-iframe] for i in tRidgeList]
outputridge=np.column_stack((tRidgeList,EridgeList,Eante,Epost))
np.savetxt('MinimizeSegmentWithRidges'+label+'_Eridge.txt', outputridge,fmt='%g %.14g %.14g %.14g', header='1)iframe 2)Eridge 3)Eante 4)Epost')

