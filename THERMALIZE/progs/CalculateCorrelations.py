#!/usr/bin/env python

import sys
import numpy as np
import argparse
from glob import glob
from math import isclose
import hoomd
from hoomd import md
import lib.module_measurements as med
import lib.module_potentials as pot
from scipy.stats import sem
import gsd.pygsd
import gsd.hoomd

parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--box_size', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='time step')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--msd', action='store_true', help='calculate mean square displacement')
parser.add_argument('--Fkt', action='store_true', help='calculate self-intermediate scattering function')
parser.add_argument('--CFF', action='store_true', help='calculate force-force correlation function')
parser.add_argument('--CFP', action='store_true', help='calculate force-momentum correlation function')
parser.add_argument('--CPP', action='store_true', help='calculate momentum-momentum correlation function')
parser.add_argument('--Cd' , action='store_true', help='calculate diagonal correlation function')
parser.add_argument('--limit_input' , type=int, required=False, default=None, help='for test runs, limits the amount of input to a small number')
args = parser.parse_args()


print(sys.argv)
print('N:',args.Natoms)
print('L:',args.box_size)
print('T:',args.temperature)
print('dt:',args.dt)
L=args.box_size
trajDIR='trajectories/'
MAXTIME=100000 # Maximum trajectory length is usually 1000, so 1e5 far more than enough


if not args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
if not args.box_size>0: raise ValueError('L ({}) must be positive'.format(args.box_size))
if not args.dt>0: raise ValueError('dt ({}) must be positive'.format(args.dt))
if not args.dt<0.1: raise ValueError('I do not believe that dt is so big ({})'.format(args.dt))
if args.temperature<0:
	raise ValueError('Temperature='+str(args.temperature)+'. Must be positive. Aborting.')
else:
	invT=1./args.temperature

def ReadAll():
	'''
	Reads all the observables from files
	'''
	def ReadTimes():
		'''
		Reads time list from the trajectory in the .gsd files
		'''
		names='./S*/{}/trajectory{}*.gsd'.format(trajDIR, args.thermostat)
		files=sorted(glob(names)) if args.limit_input==None else sorted(glob(names))[0:args.limit_input]
		ntraj=len(files)
		if 0==ntraj:
			raise NameError('No trajectory found with the name '+names)
		nt=None
		times,bad_indices=[],[]
		itraj=0
		for file in files:
			print('FILE:',file)
			with open(file, 'rb') as flow:
				HoomdFlow = gsd.pygsd.GSDFile(flow)
				hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
				nFrames = len(hoomdTraj)
				for i in range(nFrames):
					conf=hoomdTraj.read_frame(i)
					# print('iframe = ',i, '  t = ',conf.configuration.step)
					if nt==None:
						times.append(conf.configuration.step)
					elif times[i]!=conf.configuration.step:
						raise ValueError('Inconsistency in the '+str(i)+'th read time ('+str(conf.configuration.step)+' instead of '+str(times[i])+')')
			#Check nt
			if nt == None: 
				nt=nFrames
			elif nt != nFrames: 
				print('Warning: file '+file+' has nFrames='+str(nFrames)+', differently from the rest, which have nFrames='+str(nt))
				bad_indices.append(itraj)
			#Check Natoms
			if conf.particles.N != args.Natoms:
				raise ValueError('File '+file+' has and inconsistent Natoms ('+str(conf.particles.N)+' instead of '+str(args.Natoms)+')')
			itraj+=1
		if nt != len(times):
			raise ValueError('nt ('+str(nt)+') != len(times) '+str(len(times))+'')
		return ntraj,times, bad_indices

	def ReadObs(obstype):
		'''
		Reads single observable
		'''
		names='./S*/{}/{}{}*.npy'.format(trajDIR, obstype, args.thermostat)
		files=sorted(glob(names)) if args.limit_input==None else sorted(glob(names))[0:args.limit_input]
		entriesOld=None
		ntraj=len(files)
		obs, bad_indices=[], [] #List of obs, List of itw's where the trajectory is not complete
		if 0==ntraj:
			raise NameError('No trajectory found with the name '+names)
		for file in files:
			entries=0
			obstraj=[]
			fobs=open(file,'rb')
			# Loop over all the configurations of the trajectory
			for i in range(MAXTIME): 
				try: 
					conf=np.load(fobs)
					obstraj.append(conf)
					entries+=1
				except (EOFError,OSError): # I don't know why it throws also OSError when EOF is reached
					break
			fobs.close()
			# Check that Natoms is consistent
			if len(obstraj[-1]) != args.Natoms: 
				raise ValueError('File '+file+' has an inconsistent Natoms ('+str(len(obstraj[-1]))+' instead of the command-line given '+str(args.Natoms)+')')
			# Mark trajectories with inconsistent length
			if (entriesOld!=entries) and (entriesOld!=None):
				print('Warning: '+file+' has an inconsistent number of configurations ('+str(entries)+' instead of '+str(entriesOld)+')')
				bad_indices.append(len(obs))

			obs.append(obstraj)

			if entriesOld==None: 
				entriesOld=entries
		return obs, ntraj,entriesOld, bad_indices

	ntrajTimes,times,badTimes=ReadTimes()

	pos, ntrajPos, trajlenPos, badPos=ReadObs('pos')
	vel, ntrajVel, trajlenVel, badVel=ReadObs('vel')
	acc, ntrajAcc, trajlenAcc, badAcc=ReadObs('acc')
	bad_indices=set(badPos+badVel+badAcc+badTimes)
	print('incomplete indices:',bad_indices)
	for i in bad_indices:
		if i<len(pos[i]): del pos[i] # The if is because sometimes there is nothing to remove, since bad indices come from independent lists
		if i<len(vel[i]): del vel[i]
		if i<len(acc[i]): del acc[i]
	ntrajPos -= len(bad_indices)
	ntrajVel -= len(bad_indices)
	ntrajAcc -= len(bad_indices)
	print('shape(pos):',np.shape(pos))
	if 3!=np.shape(pos)[3]: raise ValueError('This should be 3D data, so make sure that the positions are 3D vectors, instead of '+str(np.shape(pos)[3])+'D')
	if 3!=np.shape(vel)[3]: raise ValueError('This should be 3D data, so make sure that the velocities are 3D vectors, instead of '+str(np.shape(vel)[3])+'D')
	if 3!=np.shape(acc)[3]: raise ValueError('This should be 3D data, so make sure that the accelarations are 3D vectors, instead of '+str(np.shape(acc)[3])+'D')
	if args.Natoms!=np.shape(pos)[2]: raise ValueError('Command-line Natoms ('+str(args.Natoms)+') does not match the one read in the Pos binary file ('+str(np.shape(pos)[2])+')')
	if args.Natoms!=np.shape(vel)[2]: raise ValueError('Command-line Natoms ('+str(args.Natoms)+') does not match the one read in the Vel binary file ('+str(np.shape(vel)[2])+')')
	if args.Natoms!=np.shape(acc)[2]: raise ValueError('Command-line Natoms ('+str(args.Natoms)+') does not match the one read in the Acc binary file ('+str(np.shape(acc)[2])+')')
	if trajlenPos!=np.shape(pos)[1]: raise ValueError('trajlenPos ('+str(trajlenPos)+') does not match the one read in the Pos binary file ('+str(np.shape(pos)[1])+')')
	if trajlenVel!=np.shape(vel)[1]: raise ValueError('trajlenVel ('+str(trajlenVel)+') does not match the one read in the Vel binary file ('+str(np.shape(vel)[1])+')')
	if trajlenAcc!=np.shape(acc)[1]: raise ValueError('trajlenAcc ('+str(trajlenAcc)+') does not match the one read in the Acc binary file ('+str(np.shape(acc)[1])+')')
	if ntrajPos!=np.shape(pos)[0]: raise ValueError('ntrajPos ('+str(ntrajPos)+') does not match the one read in the Pos binary file ('+str(np.shape(pos)[0])+')')
	if ntrajVel!=np.shape(vel)[0]: raise ValueError('ntrajVel ('+str(ntrajVel)+') does not match the one read in the Vel binary file ('+str(np.shape(vel)[0])+')')
	if ntrajAcc!=np.shape(acc)[0]: raise ValueError('ntrajAcc ('+str(ntrajAcc)+') does not match the one read in the Acc binary file ('+str(np.shape(acc)[0])+')')
	if trajlenVel  !=trajlenPos: raise ValueError('trajlenVel ('+str(trajlenVel)+') is different from trajlenPos ('+str(trajlenPos)+')')
	if trajlenVel  !=trajlenAcc: raise ValueError('trajlenVel ('+str(trajlenVel)+') is different from trajlenAcc ('+str(trajlenAcc)+')')
	if len(times)!=trajlenAcc: raise ValueError('len(times) ('+str(len(times))+') is different from trajlenAcc ('+str(trajlenAcc)+')')
	if ntrajVel  !=ntrajPos: raise ValueError('ntrajVel   ('+str(ntrajVel)  +') is different from ntrajPos ('+str(ntrajPos)+')')
	if ntrajVel  !=ntrajAcc: raise ValueError('ntrajVel   ('+str(ntrajVel)  +') is different from ntrajAcc ('+str(ntrajAcc)+')')
	if ntrajTimes!=ntrajAcc: raise ValueError('ntrajTimes ('+str(ntrajTimes)+') is different from ntrajAcc ('+str(ntrajAcc)+')')

	return (ntrajPos, times, pos, vel, acc)







def CalculateCorrelations(pos,vel,acc):
	'''
	Calculate the desired self-correlation functions
	'''
	ntw=len(pos)
	nt=len(pos[0])
	if not (ntw==len(vel   ) and ntw==len(acc   )): raise ValueError('In CalculateCorrelations, pos, vel, and acc have different ntw')
	if not (nt ==len(vel[0]) and nt ==len(acc[0])): raise ValueError('In CalculateCorrelations, pos, vel, and acc have different nt ')

	# Initialization
	if args.msd: 
		msd = np.zeros((ntw, nt), dtype=np.float64)
	if args.Fkt: 
		Fkt = np.zeros((ntw, nt), dtype=np.float64)
		n1=1; n2=3; n3=4 # Wave vector for the self-intermediate scattering function, k =[2 pi/L](n1,n2,n3) and permutations
	if args.CFF: 
		CFF = np.zeros((ntw, nt), dtype=np.float64)
	if args.CFP: 
		CFP = np.zeros((ntw, nt), dtype=np.float64)
	if args.CPP: 
		CPP = np.zeros((ntw, nt), dtype=np.float64)
	if args.Cd: 
		Cd = np.zeros((ntw, nt), dtype=np.float64)
		hoomd.context.initialize('');
		hoomd.option.set_notice_level(0)
		system = hoomd.init.read_gsd(filename='./S0/thermalized.gsd')
		if not isclose(L,system.box.Lx):
			raise ValueError('Box size L ('+str(L)+') is inconsistent with the read sample configuration ('+str(system.box.Lx)+')')
		pair=pot.LJ(md.nlist.cell(), type='KA' if args.Natoms>500 else 'KAshort')
		snapA=system.take_snapshot()
		snapB=system.take_snapshot()

	# Calculate the correlation functions
	obs={}
	for itw in range(ntw):
		print('itw:',itw)
		initialPositions=np.copy(pos[itw][0]) #deep copy actually not needed bc there arrays are not modified
		if args.msd:
			msd[itw] = np.array([med.PeriodicSquareDistance(pos[itw][it], initialPositions, L) for it in range(nt)])/args.Natoms
		if args.Fkt:
			Fkt[itw] = np.array([med.ComputeFkt(n1, n2, n3, L, med.PeriodicDisplacement(pos[itw][it], initialPositions, L)) for it in range(nt)])
		if args.CFF:
			CFF[itw] = np.array([np.mean([np.inner(acc[itw][0][atom],acc[itw][it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
		if args.CFP:
			CFP[itw] = np.array([np.mean([np.inner(acc[itw][0][atom],vel[itw][it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
		if args.CPP:
			CPP[itw] = np.array([np.mean([np.inner(vel[itw][0][atom],vel[itw][it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
		if args.Cd :
			snapA.particles.position[:] = initialPositions[:] #Deep copy not needed bc there arrays are not modified
			for it in range(nt):
				print('\rit:',it)
				snapB.particles.position[:] = pos[itw][it]
				Cd [itw][it]=pair.Cd(snapA=snapA, snapB=snapB, beta=invT)

	# Save correlation functions on an array
	if args.msd: obs['msd']={'mean': np.mean(msd,axis=0), 'err': sem(msd, axis=0)}
	if args.Fkt: obs['Fkt']={'mean': np.mean(Fkt,axis=0), 'err': sem(Fkt, axis=0)}
	if args.CFF: obs['CFF']={'mean': np.mean(CFF,axis=0), 'err': sem(CFF, axis=0)}
	if args.CFP: obs['CFP']={'mean': np.mean(CFP,axis=0), 'err': sem(CFP, axis=0)}
	if args.CPP: obs['CPP']={'mean': np.mean(CPP,axis=0), 'err': sem(CPP, axis=0)}
	if args.Cd : obs['Cd' ]={'mean': np.mean(Cd, axis=0), 'err': sem(Cd , axis=0)}

	return obs


def WriteCorrelations(times):
	'''
	Writes to file the measured correlations, both in binary (for further analysis) and in text (for plots)
	'''
	times=np.array(times)*args.dt
	np.save('times_{}.npy'.format(args.thermostat), times)
	if args.msd:
		np.savetxt('msd_{}.txt'.format(args.thermostat), np.column_stack((times, obs['msd']['mean'], obs['msd']['err'])), fmt=['%.14g','%.14g','%.14g'], header='time msd err')
		np.save('msd_{}.npy'.format(args.thermostat), obs['msd'])
	if args.Fkt:
		np.savetxt('Fkt_{}.txt'.format(args.thermostat), np.column_stack((times, obs['Fkt']['mean'], obs['Fkt']['err'])), fmt=['%.14g','%.14g','%.14g'], header='time Fkt err')
		np.save('Fkt_{}.npy'.format(args.thermostat), obs['Fkt'])
	if args.CFF:
		np.savetxt('CFF_{}.txt'.format(args.thermostat), np.column_stack((times, obs['CFF']['mean'], obs['CFF']['err'])), fmt=['%.14g','%.14g','%.14g'], header='time CFF err')
		np.save('CFF_{}.npy'.format(args.thermostat), obs['CFF'])
	if args.CFP:
		np.savetxt('CFP_{}.txt'.format(args.thermostat), np.column_stack((times, obs['CFP']['mean'], obs['CFP']['err'])), fmt=['%.14g','%.14g','%.14g'], header='time CFP err')
		np.save('CFP_{}.npy'.format(args.thermostat), obs['CFP'])
	if args.CPP:
		np.savetxt('CPP_{}.txt'.format(args.thermostat), np.column_stack((times, obs['CPP']['mean'], obs['CPP']['err'])), fmt=['%.14g','%.14g','%.14g'], header='time CPP err')
		np.save('CPP_{}.npy'.format(args.thermostat), obs['CPP'])
	if args.Cd :
		np.savetxt('Cd_{}.txt'.format(args.thermostat) , np.column_stack((times, obs['Cd' ]['mean'], obs['Cd' ]['err'])), fmt=['%.14g','%.14g','%.14g'], header='time Cd err')
		np.save('Cd_{}.npy'.format(args.thermostat), obs['Cd'])
	return




(ntw,times,pos,vel,acc)=ReadAll()
obs=CalculateCorrelations(pos,vel,acc)
WriteCorrelations(times)



