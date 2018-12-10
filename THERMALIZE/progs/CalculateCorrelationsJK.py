#!/usr/bin/env python
#
# This program does the same as CalculateCorrelations (i.e. it calculates the trivial
# correlation functions: msd, Fkt, Cd, CFF, CFP, CPP), but it separates it in jackknife blocks.
# This has the advantage of allowing error computing and of using less memory (since blocks
# are made during reding)


import sys,os
import numpy as np
import argparse
from glob import glob
# from math import isclose
import hoomd
from hoomd import md
import lib.module_measurements as med
# import lib.module_potentials as pot
# from scipy.stats import sem
import gsd.pygsd
import gsd.hoomd
from matplotlib import pyplot as plt


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
parser.add_argument('--ntCd', type=int, required=False, default=30, help='Number of times for the calculation of Cd')
parser.add_argument('--limit_input' , type=int, required=False, default=None, help='for test runs, limits the amount of input to a small number')
parser.add_argument('--lblo', type=int, required=False, default='10', help='Length of jackknife blocks')
args = parser.parse_args()


print(sys.argv)
print('N:',args.Natoms)
print('L:',args.box_size)
print('T:',args.temperature)
print('dt:',args.dt)
L=args.box_size
trajDIR='trajectories/'
maxnsam=100
maxntraj=1000

readPos = True if (args.msd or args.Fkt or args.Cd) else False
readVel = True if (args.CPP or args.CFP)            else False
readAcc = True if (args.CFF or args.CFP)            else False

if not args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
if not args.box_size>0: raise ValueError('L ({}) must be positive'.format(args.box_size))
if not args.dt>0: raise ValueError('dt ({}) must be positive'.format(args.dt))
if not args.dt<0.1: raise ValueError('I do not believe that dt is so big ({})'.format(args.dt))
if args.temperature<0:
	raise ValueError('Temperature='+str(args.temperature)+'. Must be positive. Aborting.')
else:
	invT=1./args.temperature


msd, Fkt, CFF, CFP, CPP = {'blocks': []}, {'blocks': []}, {'blocks': []}, {'blocks': []}, {'blocks': []}

def ReadAllBlocks():

	count=0
	iblo=-1
	for isam in range(maxnsam):
		samdir="./S{}/{}".format(isam,trajDIR)
		if not os.path.isdir(samdir): 
			# print("Skipping sample {}".format(isam))
			continue
		for itraj in range(maxntraj):
			namegsd="{}/trajectory{}{}.gsd".format(samdir, args.thermostat, itraj)
			namePos="{}/pos{}{}.npy".format(samdir, args.thermostat, itraj)
			nameVel="{}/vel{}{}.npy".format(samdir, args.thermostat, itraj)
			nameAcc="{}/acc{}{}.npy".format(samdir, args.thermostat, itraj)
			if not os.path.isfile(namegsd):
				# print('Skipping sample {} traj {} because {} does not exist'.format(isam, itraj, namegsd))
				continue
			if not os.path.isfile(namePos):
				# print('Skipping sample {} traj {} because {} does not exist'.format(isam, itraj, namePos))
				continue
			if not os.path.isfile(nameVel):
				# print('Skipping sample {} traj {} because {} does not exist'.format(isam, itraj, nameVel))
				continue
			if not os.path.isfile(nameAcc):
				# print('Skipping sample {} traj {} because {} does not exist'.format(isam, itraj, nameAcc))
				continue
			print('\risam = ',isam, ' itraj = ',itraj, end='')

			if not (os.path.isfile(namegsd) and os.path.isfile(namePos) and os.path.isfile(nameVel) and os.path.isfile(nameAcc)):
				print('Skipping sample {} traj {}'.format(isam, itraj))
				continue

			#Reat Times
			if count==0:
				times=ReadTimesSample(namegsd)
				nt=len(times)
			else:
				if not np.array_equal(ReadTimesSample(namegsd),times):
					print('Inconsistent times - skipping...')
					sys.exit('fine prova')
					continue

			if count%args.lblo==0:
				iblo+=1
				print('iblo: ',iblo)
				n_in_blo=0
				if args.msd: msd_blo=np.zeros(nt, dtype=np.float64)
				if args.Fkt: Fkt_blo=np.zeros(nt, dtype=np.float64)
				if args.CFF: CFF_blo=np.zeros(nt, dtype=np.float64)
				if args.CFP: CFP_blo=np.zeros(nt, dtype=np.float64)
				if args.CPP: CPP_blo=np.zeros(nt, dtype=np.float64)

			if readPos:
				pos=ReadTraj(namePos, nt)
			if readVel:
				vel=ReadTraj(nameVel, nt)
			if readAcc:
				acc=ReadTraj(nameAcc, nt)

			if args.msd:
				this_msd = np.array([med.PeriodicSquareDistance(pos[it], pos[0], L) for it in range(nt)])/args.Natoms
				msd_blo+=this_msd
			if args.Fkt: 
				n1=1; n2=3; n3=4 # Wave vector for the self-intermediate scattering function, k =[2 pi/L](n1,n2,n3) and permutations
				this_Fkt = np.array([med.ComputeFkt(n1, n2, n3, L, med.PeriodicDisplacement(pos[it], pos[0], L)) for it in range(nt)])
				Fkt_blo+=this_Fkt
			if args.CFF: 
				this_CFF = np.array([np.mean([np.inner(acc[0][atom],acc[it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
				CFF_blo+=this_CFF
			if args.CFP: 
				this_CFP = np.array([np.mean([np.inner(acc[0][atom],vel[it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
				CFP_blo+=this_CFP
			if args.CPP: 
				this_CPP = np.array([np.mean([np.inner(vel[0][atom],vel[it][atom]) for atom in range(args.Natoms)]) for it in range(nt)])/3.
				CPP_blo+=this_CPP
			if args.Cd: 
				raise NotImplementedError('ReadAllBlocks(): I have not yet implemented the calculation of Cd in this program')

			n_in_blo+=1
			count+=1

			if count%args.lblo==0:
				assert(n_in_blo<=args.lblo)
				if args.msd: msd['blocks'].append(msd_blo/n_in_blo)
				if args.Fkt: Fkt['blocks'].append(Fkt_blo/n_in_blo)
				if args.CFF: CFF['blocks'].append(CFF_blo/n_in_blo)
				if args.CFP: CFP['blocks'].append(CFP_blo/n_in_blo)
				if args.CPP: CPP['blocks'].append(CPP_blo/n_in_blo)

	print('')

	#The last, incomplete block
	if count%args.lblo!=0:
		if args.msd: msd['blocks'].append(msd_blo/n_in_blo)
		if args.Fkt: Fkt['blocks'].append(Fkt_blo/n_in_blo)
		if args.CFF: CFF['blocks'].append(CFF_blo/n_in_blo)
		if args.CFP: CFP['blocks'].append(CFP_blo/n_in_blo)
		if args.CPP: CPP['blocks'].append(CPP_blo/n_in_blo)

	if args.msd: 
		nblo=len(msd['blocks']); nm1=nblo-1
		msd['blocks']  = np.array(msd['blocks'])
		msd['mean'  ]  = msd['blocks'].mean(axis=0)
		msd['blocksJK']= np.array([msd['mean']*nblo - block for block in msd['blocks']])/nm1
		msd['errJK']   = np.sqrt(nm1*(np.square(msd['blocksJK']).mean(axis=0) - msd['mean']*msd['mean']))
	if args.Fkt: 
		nblo=len(Fkt['blocks']); nm1=nblo-1
		Fkt['blocks']=np.array(Fkt['blocks'])
		Fkt['mean'  ]=Fkt['blocks'].mean(axis=0)
		Fkt['blocksJK']= np.array([Fkt['mean']*nblo - block for block in Fkt['blocks']])/(nblo-1)
		Fkt['errJK']   = np.sqrt(nm1*(np.square(Fkt['blocksJK']).mean(axis=0) - Fkt['mean']*Fkt['mean']))
	if args.CFF: 
		nblo=len(CFF['blocks']); nm1=nblo-1
		CFF['blocks']=np.array(CFF['blocks'])
		CFF['mean'  ]=CFF['blocks'].mean(axis=0)
		CFF['blocksJK']= np.array([CFF['mean']*nblo - block for block in CFF['blocks']])/nm1
		CFF['errJK']   = np.sqrt(nm1*(np.square(CFF['blocksJK']).mean(axis=0) - CFF['mean']*CFF['mean']))
	if args.CFP: 
		nblo=len(CFP['blocks']); nm1=nblo-1
		CFP['blocks']=np.array(CFP['blocks'])
		CFP['mean'  ]=CFP['blocks'].mean(axis=0)
		CFP['blocksJK']= np.array([CFP['mean']*nblo - block for block in CFP['blocks']])/nm1
		CFP['errJK']   = np.sqrt(nm1*(np.square(CFP['blocksJK']).mean(axis=0) - CFP['mean']*CFP['mean']))
	if args.CPP: 
		nblo=len(CPP['blocks']); nm1=nblo-1
		CPP['blocks']=np.array(CPP['blocks'])
		CPP['mean'  ]=CPP['blocks'].mean(axis=0)
		CPP['blocksJK']= np.array([CPP['mean']*nblo - block for block in CPP['blocks']])/nm1
		CPP['errJK']   = np.sqrt(nm1*(np.square(CPP['blocksJK']).mean(axis=0) - CPP['mean']*CPP['mean']))


	return times

def ReadTraj(filename, ntimes):
	obstraj=[]
	fobs=open(filename,'rb')
	# Loop over all the configurations of the trajectory
	for i in range(ntimes): 
		try: 
			conf=np.load(fobs)
			obstraj.append(conf)
		except (EOFError,OSError): # I don't know why it throws also OSError when EOF is reached
			raise IndexError('ReadObs is reading too many lines')
			break
	fobs.close()
	return np.array(obstraj)


def ReadTimesSample(filename):
	times=[]
	with open(filename, 'rb') as flow:
		HoomdFlow = gsd.pygsd.GSDFile(flow)
		hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
		nFrames = len(hoomdTraj)
		for i in range(nFrames):
			conf=hoomdTraj.read_frame(i)
			# print('iframe = ',i, '  t = ',conf.configuration.step)
			times.append(conf.configuration.step)
			if conf.particles.N != args.Natoms:
				raise ValueError('File '+file+' has and inconsistent Natoms ('+str(conf.particles.N)+' instead of '+str(args.Natoms)+')')
	return np.array(times)





def WriteCorrelations(times):
	'''
	Writes to file the measured correlations, both in binary (for further analysis) and in text (for plots)
	'''
	times=np.array(times)*args.dt
	np.save('times_{}.npy'.format(args.thermostat), times)
	if args.msd:
		np.savetxt('msdJK_{}.txt'.format(args.thermostat), np.column_stack((times, msd['mean'], msd['errJK'])), fmt=['%.14g','%.14g','%.14g'], header='time msd err')
		np.save('msdJK_{}.npy'.format(args.thermostat), msd)
	if args.Fkt:
		np.savetxt('FktJK_{}.txt'.format(args.thermostat), np.column_stack((times, Fkt['mean'], Fkt['errJK'])), fmt=['%.14g','%.14g','%.14g'], header='time Fkt err')
		np.save('FktJK_{}.npy'.format(args.thermostat), Fkt)
	if args.CFF:
		np.savetxt('CFFJK_{}.txt'.format(args.thermostat), np.column_stack((times, CFF['mean'], CFF['errJK'])), fmt=['%.14g','%.14g','%.14g'], header='time CFF err')
		np.save('CFFJK_{}.npy'.format(args.thermostat), CFF)
	if args.CFP:
		np.savetxt('CFPJK_{}.txt'.format(args.thermostat), np.column_stack((times, CFP['mean'], CFP['errJK'])), fmt=['%.14g','%.14g','%.14g'], header='time CFP err')
		np.save('CFPJK_{}.npy'.format(args.thermostat), CFP)
	if args.CPP:
		np.savetxt('CPPJK_{}.txt'.format(args.thermostat), np.column_stack((times, CPP['mean'], CPP['errJK'])), fmt=['%.14g','%.14g','%.14g'], header='time CPP err')
		np.save('CPPJK_{}.npy'.format(args.thermostat), CPP)
	return


times=ReadAllBlocks()
WriteCorrelations(times)




