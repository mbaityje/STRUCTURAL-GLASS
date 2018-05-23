#!/usr/bin/env python
################################################################
#
# DESCRIPTION
# This example reads a trajectory and finds the related
# trajectory of the inherent structures.
# Also, the ridge energies between ISs are recorded.
#
################################################################

#MEMO:
#-Le posizioni devono essere in doppia precisione
#-Bisogna trovare il modo di fare accordare chunks successivi
#

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import argparse
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import gsd.pygsd
import gsd.hoomd
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import module_measurements as med
import os.path

print("#++++++++++++++++#")
print("#----------------#")
print("# BisectChunk.Py #")
print("#----------------#")
print("#++++++++++++++++#")

################################################################
#
# READ ARGUMENTS
# 
################################################################

#Integrator parameters
alphaFIRE=0.99
ftolFIRE=1e-5
EtolFIRE=1e-10
wtolFIRE=1e-5
minstepsFIRE=100


#Start hoomd
print("Initialize hoomd context\n")
generalContext=hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd
#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', help='name of the trajectory file .gsd')
parser.add_argument('--ichunk', type=int, required=True, help='index of the chunk')
parser.add_argument('--tchunk', type=int, required=True, help='number of MD steps in this chunk')
parser.add_argument('-l','--label', required=False, default='', help='label for distinguishing runs and continuations')
parser.add_argument('--deltaE', type=float, required=False, default=1e-6, help='energy difference in order to consider different two configurations')
parser.add_argument('--skiprows', type=int, required=False, default=0, help='how many rows we can skip when reading elist (we need only the last one)')
parser.add_argument('--doridge', action='store_true', help='If activated, calculate ridge between ISs.')
parser.add_argument('--moreobs', action='store_true', help='If activated, calculate nmoved, msd and q.')
parser.add_argument('--verbose', action='store_true', help='If activated, print extra information to stdout.')
args = parser.parse_args(more_arguments)


filename=args.filename
del parser

print("filename = ",args.filename)
print("ichunk   = ",args.ichunk)
print("tchunk   = ",args.tchunk)
print("deltaE   = ",args.deltaE)
print("doRidge  = ",args.doridge)
print("moreobs  = ",args.moreobs)
print("verbose  = ",args.verbose)



################################################################
#
# READ TRAJECTORY
# 
################################################################

with open(args.filename, 'rb') as flow:
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
	
	Nframes = len(hoomdTraj)
	print('There are Nframes=', Nframes, 'in the file.')
	if Nframes != args.tchunk: raise ValueError('Nframes != args.tchunk')

	posizioni=[hoomdTraj[i].particles.position[:] for i in range(Nframes)]
	HoomdFlow.close()


################################################################
# 
# Initialize
#
################################################################
system = hoomd.init.read_gsd(filename=args.filename)
t0=hoomd.get_step()
print("Initial time:", t0)
snap_ini=system.take_snapshot()
snap_final=system.take_snapshot()
snap_final.particles.position[:]=posizioni[Nframes-1]

#If it's the first chunk, there is no list of energies.
#Otherwise, we open it and make sure that the time step is consistent.
if(args.ichunk>0):
	elist_old=np.loadtxt('elistIS.txt',skiprows=args.skiprows)
	if int(elist_old[len(elist_old)-1][0])!=t0-1: raise ValueError('Discrepacy with the previous chunk:\n elist_old[len(elist_old)-1][0])='+str(int(elist_old[len(elist_old)-1][0]))+', t0-1='+str(t0-1))
	del elist_old

################################################################
# 
# Set potential and analyzer
#
################################################################
NeighborsList = md.nlist.cell()
potential=pot.LJ(NeighborsList,type="KAshort") #myLJpair is now an attribute of potential. To call it: potential.GetLJpair()
analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
analyzer = hoomd.analyze.log(filename=None, quantities=analyzer_quantities, period=1)

################################################################
# 
# Function declaration
#
################################################################

# This function is commented because I want to avoid using it, but it can be useful for debugging
# def PotEn(mode,analyzer,dt=0.0025):
# 	'''
# 	Returns system's energy if the analyzer is set
# 	'''
# 	mode.dt=1e-24
# 	hoomd.run(2)
# 	U=analyzer.query('potential_energy')
# 	mode.dt=dt
# 	return U

def Minimize(snap):
	system.restore_snapshot(snap)
	fire.cpp_integrator.reset()
	if not integrator_fire.enabled: integrator_fire.enable()
	while not(fire.has_converged()):
		hoomd.run(100)
	eIS=analyzer.query('potential_energy')
	# snapnew=system.take_snapshot()
	integrator_fire.disable()
	return eIS


def Bisect(t1, snap1, t2, snap2, e1=None, doRidge=False):
	assert(t2>t1)
	# Energies of the inherent structures
	# If e1 was given as an argument, no need to caclulate it again
	if (e1==None):
		e1=Minimize(snap1)
	e2=Minimize(snap2)
	ediff=np.abs(e1-e2)
	print("Le due eIS:",t1,e1," ,   ",t2,e2,"       diff=",e1-e2)

	#If they are next to each other, I record them no matter their energy
	if (t2-t1)==1:
		if doRidge and (ediff>args.deltaE): 
			eRidge,snapRidge=CalculateRidge(snap1, snap2, e1, e2, L, verbose=False, time=t0+t2)
			if eRidge != None:
				ridgetime=t0+t2-0.5
				eRidgelist[ridgetime] = eRidge
				if args.moreobs:
					nmovedlist['12'][ridgetime] = med.HowManyMovedPos(    snap1.particles.position,     snap2.particles.position, L)
					nmovedlist['1r'][ridgetime] = med.HowManyMovedPos(    snap1.particles.position, snapRidge.particles.position, L)
					nmovedlist['r2'][ridgetime] = med.HowManyMovedPos(snapRidge.particles.position,     snap2.particles.position, L)
					msdlist['12'][ridgetime] = med.PeriodicSquareDistance(    snap1.particles.position,     snap2.particles.position, L)
					msdlist['1r'][ridgetime] = med.PeriodicSquareDistance(    snap1.particles.position, snapRidge.particles.position, L)
					msdlist['r2'][ridgetime] = med.PeriodicSquareDistance(snapRidge.particles.position,     snap2.particles.position, L)
					qlist['12'][ridgetime] = med.OverlapPos(    snap1.particles.position,     snap2.particles.position, L)
					qlist['1r'][ridgetime] = med.OverlapPos(    snap1.particles.position, snapRidge.particles.position, L)
					qlist['r2'][ridgetime] = med.OverlapPos(snapRidge.particles.position,     snap2.particles.position, L)
		elist[t0+t2]=e2

		return

	#If they are the same I save and finish
	if(ediff<args.deltaE):
		elist[t0+t2]=e2
		return
	# If they are different, I find the intermediate time (that is necessarily
	# different because I made sure a few lines above), and I bisect both the first
	# and the second interval.
	else:
		t12=int(0.5*(t1+t2))
		snap12=system.take_snapshot()
		snap12.particles.position[:]=posizioni[t12][:]
		Bisect(t1,snap1,t12,snap12,e1, doRidge=doRidge)
		Bisect(t12,snap12,t2,snap2, doRidge=doRidge)


def ConfBisect(snap1, snap2, eis1, eis2, L, dmax=0.002, verbose=False):
	"""
	This function takes two thermal snapshots that end in two separate IS,
	and makes a linear interpolation between the positions of the two.
	Through bisection, the point that separates the two basins is found.
	#0.002 is about half the typical distance between confs at subsequent time steps w/ dt=0.0025
	"""
	if verbose: print('···ConfBisect···')
	if verbose: print('eis1 = ',eis1,'; eis2 = ',eis2)

	if np.abs(eis1-eis2)<args.deltaE: raise ValueError('ConfBisect ERROR: the two starting ISs the same one [abs(eis1-eis2)<args.deltaE]')
	dstart=0.5*dmax
	Natoms=snap1.particles.N
	dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms #the box is cubic
	snap12=LinearConfInterpolation(snap1, snap2, L)

	eis12=Minimize(snap12)
	if verbose: 
		print('snap1 .particles.position[0] = ',snap1.particles.position[0])
		print('snap2 .particles.position[0] = ',snap2.particles.position[0])
		print('snap12.particles.position[0] = ',snap12.particles.position[0])


	count=0
	maxcount=100
	while dist12>dstart:
		if verbose: print('eis1 = ',eis1,'; eis2 = ',eis2,'; eis12 = ',eis12)
		if np.abs(eis1-eis12) <= args.deltaE: #If snap12 belongs to snap1, snap1=snap12
			snap1.particles.position[:]=snap12.particles.position[:]
			eis1=eis12
		elif np.abs(eis2-eis12) <= args.deltaE: #If snap12 belongs to snap2, snap2=snap12
			snap2.particles.position[:]=snap12.particles.position[:]
			eis2=eis12
		else: #If snap12 does not belong to either, we throw a warning and change snap2
			print("ConfBisect: found an intermediate IS while searching the TS (Eis12=%.14g)"%eis12)
			snap2.particles.position[:]=snap12.particles.position[:]
			eis2=eis12
		dist12=med.PeriodicDistance(snap1.particles.position, snap2.particles.position, L).sum()/Natoms
		if verbose: print("dist12=",dist12)
			
		count+=1
		if  np.abs(eis1-eis2)<args.deltaE: raise ValueError('ConfBisect ERROR: the two ISs became the same one [abs(eis1-eis2)<args.deltaE]')
		snap12=LinearConfInterpolation(snap1, snap2, L)
		if verbose: print('snap12.particles.position[0] = ',snap12.particles.position[0])
		eis12=Minimize(snap12)

		if count>maxcount:
			raise RecursionError("ConfBisect ERROR: the interpolation bisection is not converging.")
	return snap1,snap2,snap12,eis1,eis2,eis12,dist12

def LinearConfInterpolation(snap1, snap2, box_size):
	'''
	returns a linear interpolation between snap1 and snap2
	new_positions = (snap1 + snap2)/2
	'''
	snap12=system.take_snapshot()
	snap12.particles.position[:] = med.PeriodicIntermPoints(
		np.array(snap1.particles.position,dtype=np.float64),
		np.array(snap2.particles.position,dtype=np.float64),
		box_size)
	return snap12

def CalcF():
	'''Calculates the forces of the system through the integrator
	Mind that if the dynamics were not run, this gives zero'''
	F=np.array([system.particles[i].net_force for i in range(Natoms)])
	return F/Natoms

def CalcAcc(snapshot):
	'''Returns the accelerations of the snapshot
	Mind that some snapshots don't have information on the accelerations'''
	acc=np.array(snapshot.particles.acceleration,dtype=np.float64)
	return acc/Natoms

def SetupAnalyzer(logname=None, period=None, quantities=['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum']):
	analyzer = hoomd.analyze.log(filename=logname, \
								 quantities=quantities, period=period, \
								 header_prefix = '#', \
								 overwrite=False,
								 phase=0)
	return analyzer

def CalculateRidge(snapT1, snapT2, Eis1, Eis2, L, verbose=False, dtFIRE=0.0025, time=None):
	''' Calculates energy at the ridge between two snapshots '''
	if verbose: print("--- Calculate Ridge ---\n- Eis1:",Eis1,' Eis2:',Eis2)

	dmax=0.0001 #0.004 is about the typical distance between confs at subsequent time steps with dt=0.0025
	nsteps=1
	snapT1.particles.velocity[:]=np.zeros((Natoms, 3))
	snapT2.particles.velocity[:]=np.zeros((Natoms, 3))

	#snapis are not inherent structures, but the gradually approach them
	snapis1,snapis2,snapis12,Eis1,Eis2,Eis12,dist12=ConfBisect(snapT1, snapT2, Eis1, Eis2, L, dmax=dmax, verbose=False)

	if verbose: print("- Eis1=",Eis1,'; Eis2=',Eis2,' (after ConfBisect)')

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
	# if verbose: print("+  eis1 = ",Minimize(snapis1),"\teis2 = ",Minimize(snapis2))

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

		if verbose: print("iter: ",iter," dist12=",dist12,"eis1: %.14f"%eis1," eis2: %.14f"%eis2)

		if dist12>dmax:
			'''
			With the linear interpolations the found barrier is occasionally lower than one of the ISs, because of nonlinearities in the path.
			# snapRidge=LinearConfInterpolation(snapis1_old, snapis2_old, L)
			# Eridge=potential.CalculateEnergySlower(snapRidge)
			For this reason, I just take the energy of one of the two confs at the previous step.
			'''
			if eis1>eis2: 
				Eridge=ConsistentRidge(eis1_old,Eis1,Eis2,args.deltaE,EtolFIRE)
				snapRidge=snapis1_old
			else : 
				Eridge=ConsistentRidge(eis2_old,Eis1,Eis2,args.deltaE,EtolFIRE)
				snapRidge=snapis2_old
			print("Eridge = ",Eridge)
			break

		snapis1_old=snapis1
		snapis2_old=snapis2
		eis1_old=eis1
		eis2_old=eis2

	if iter==niter-1: 
		ridgeLog.Exception(time)
		Eridge=None
		snapRidge=None
	return Eridge,snapRidge


class RidgeConvergence:
	def __init__(self, niter, ichunk, tolerant=None, label=None):
		self.tolerant = False if tolerant==None else tolerant
		self.niter=niter
		self.label = '' if label==None else label
		self.ichunk=ichunk
		if self.tolerant==True:
			self.CreateLog()
		return

	def Exception(self, time=None):
		if self.tolerant:
			print('CalculateRidge did not converge after '+str(self.niter)+' steps')
			self.UpdateLog(time)
		else:
			raise RecursionError('CalculateRidge did not converge after '+str(self.niter)+' steps')
		return

	def CreateLog(self):
		self.filename='ridge'+self.label+'.log'
		mode = 'w' if self.ichunk==0 else 'a'
		out = open(self.filename   , mode)
		if self.ichunk==0:
			print('This log tracks the times the ridge was not found after ',self.niter,' iterations (NC: Non Converged)', file=out)
		print('ichunk = ',self.ichunk)
		out.close()
		return

	def UpdateLog(self, time):
		out = open(self.filename , 'a')
		print('NC',file=out, end=' ')
		if time==None:	print(''  , file=out)
		else:			print(time, file=out)
		out.close()
		return


def ConsistentRidge(Eridge, Eis1, Eis2, thres, EtolFIRE=0):
	'''
	If Eridge>Eis1 and Eridge>Eis2 everything is consistent.
	For precision reasons, it might happen that we must use the condition Eridge+thres>Eis.
	If neither this condition is met, something is wrong.
	'''
	if Eridge>Eis1:
		if Eridge>Eis2: 
			return Eridge
		elif Eridge+thres+EtolFIRE>Eis2:
			return Eridge+thres+EtolFIRE
		else:
			sys.exit('Inconsistent Eridge<Eis.\nEridge='+str(Eridge)+'\tEis2='+str(Eis2))
	elif Eridge+thres+EtolFIRE>Eis1:
		return Eridge+thres+EtolFIRE
	else:
		sys.exit('Inconsistent Eridge<Eis.\nEridge='+str(Eridge)+'\tEis1='+str(Eis1))

################################################################
# 
# Bisection
#
################################################################
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=alphaFIRE, ftol=ftolFIRE, Etol=EtolFIRE, wtol=wtolFIRE, min_steps=minstepsFIRE)
integrator_fire = md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)

#List of energies
elist={t0:Minimize(snap_ini)}
if args.doridge:
	eRidgelist={} # List of ridge energies
	if args.moreobs:
		nmovedlist ={'12':{}, '1r':{}, 'r2':{}} # How many particles moved around a ridge (12: the two IS, 1r: initial IS and ridge, r2: ridge and final IS)
		msdlist    ={'12':{}, '1r':{}, 'r2':{}} # Mean square displacement around a ridge (12: the two IS, 1r: initial IS and ridge, r2: ridge and final IS)
		qlist={'12':{}, '1r':{}, 'r2':{}} # Overlap around a ridge (12: the two IS, 1r: initial IS and ridge, r2: ridge and final IS)
	niter=10000
	ridgeLog=RidgeConvergence(niter,tolerant=True, ichunk=args.ichunk, label=args.label)



Bisect(0,snap_ini,Nframes-1,snap_final, e1=None, doRidge=args.doridge)


################################################################
# 
# Cleanup
#
################################################################
if integrator_fire.enabled: integrator_fire.disable()

################################################################
# 
# Output
#
################################################################
mode='wb' if args.ichunk==0 else 'ab'

#list of inherent structure energies
outIS    = open('elistIS'+args.label+'.txt', mode)
headerIS = 'time Eis' if args.ichunk==0 else ''
outputIS = np.column_stack( ( list(elist.keys()), list(elist.values() )) )
np.savetxt(outIS, outputIS,fmt='%d %.14g', header=headerIS)
outIS.close()

#list of ridge energies
if args.doridge:
	outRidge    = open('elistRidge'+args.label+'.txt', mode)
	Eante=[elist[int(t)]   for t in list(eRidgelist.keys())]
	Epost=[elist[int(t)+1] for t in list(eRidgelist.keys())]
	if args.moreobs:
		headerRidge='time Eridge Eante Epost n12 n1r nr2 msd12 msd1r msdr2 q12 q1r qr2' if args.ichunk==0 else ''
		np.savetxt(outRidge, np.column_stack(
			( 	list(eRidgelist.keys()), 
			list(eRidgelist.values()), 
			Eante, 
			Epost,
			list(nmovedlist['12'].values()),
			list(nmovedlist['1r'].values()),
			list(nmovedlist['r2'].values()),
			list(msdlist['12'].values()),
			list(msdlist['1r'].values()),
			list(msdlist['r2'].values()),
			list(qlist['12'].values()),
			list(qlist['1r'].values()),
			list(qlist['r2'].values())	)),
			fmt='%.1f %.14g %.14g %.14g %d %d %d %.14g %.14g %.14g %.14g %.14g %.14g', header=headerRidge)
	else:
		headerRidge='time Eridge Eante Epost' if args.ichunk==0 else ''
		np.savetxt(outRidge, np.column_stack(
			( 	list(eRidgelist.keys()), 
			list(eRidgelist.values()), 
			Eante, 
			Epost )),
			fmt='%.1f %.14g %.14g %.14g', header=headerRidge)
		outRidge.close()
