#!/usr/bin/env python

import sys
import numpy as np
import argparse, time, itertools
from glob import glob
from matplotlib import pyplot as plt
from math import isclose
import hoomd
from hoomd import md
import lib.module_measurements as med
import lib.module_potentials as pot
from scipy.stats import sem
from scipy import integrate
import gsd.pygsd
import gsd.hoomd

parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--box_size', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('-n', type=int, nargs=3, required=False, default=[1,3,4], help='[n1,n2,n3]: k=(2pi/L)[n1,n2,n3]')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--ntCd', type=int, required=False, default=25, help='Number of times for the calculation of correlation')
parser.add_argument('--limit_input' , type=int, required=False, default=None, help='for test runs, limits the amount of input to a small number')
args = parser.parse_args()

L=args.box_size
args.n=np.array(args.n, dtype=int)

if 1: #All permutations of k components
	sqn=np.inner(args.n,args.n)
	all_ni=np.concatenate((args.n,-args.n))
	allperm=list(itertools.permutations(all_ni,3))
	all_n=np.unique(np.array(list(filter(lambda x: np.inner(x,x)==sqn, allperm))), axis=0 )
else: #Only some permutations of k components
	all_n=np.unique(list(itertools.permutations(args.n)), axis=0)
all_k=np.float64(2*np.pi/L)*all_n

print(sys.argv)
print('N:',args.Natoms)
print('L:',args.box_size)
print('T:',args.temperature)
print('n:',args.n)
print('nk:',len(all_n))



trajDIR='trajectories/'
MAXTIME=10000 # Maximum trajectory length is usually 1e3

readPos = True
readVel = False
readAcc = True

if not args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
if not args.box_size>0: raise ValueError('L ({}) must be positive'.format(args.box_size))
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
		names='../S*/{}/trajectory{}[0-9][0-9]*.gsd'.format(trajDIR, args.thermostat) #Exclude the first 10 trajectories - just a check
		all_names=sorted(glob(names))
		np.random.seed(seed=12345)
		np.random.shuffle(all_names)
		# raise KeyboardInterrupt
		files=all_names if args.limit_input==None else all_names[0:args.limit_input]
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
		names='../S*/{}/{}{}*.npy'.format(trajDIR, obstype, args.thermostat)
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

	def checks(xxx, ntraj, trajlen, name):
		''' xxx=pos,vel,acc'''
		print('shape({}): {}'.format(name, np.shape(xxx)))
		if 3!=np.shape(xxx)[3]: raise ValueError('This should be 3D data, so make sure that the {name} are 3D vectors, instead of {np.shape(xxx)[3]}D')
		if args.Natoms!=np.shape(xxx)[2]: raise ValueError('Command-line Natoms ('+str(args.Natoms)+') does not match the one read in the {name} binary file ({str(np.shape(xxx)[2])})')
		if trajlen!=np.shape(xxx)[1]: raise ValueError('trajlen ('+str(trajlen)+') does not match the one read in the {name} binary file ('+str(np.shape(xxx)[1])+')')
		if ntraj!=np.shape(xxx)[0]: raise ValueError('ntraj ('+str(ntraj)+') does not match the one read in the {name} binary file ('+str(np.shape(xxx)[0])+')')
		if len(times)  !=trajlen: raise ValueError('len(times) ('+str(len(times))+') is different from trajlen{name} ('+str(trajlen)+')')
		if ntrajTimes!=ntraj: raise ValueError('ntrajTimes ('+str(ntrajTimes)+') is different from ntraj{name} ('+str(ntraj)+')')


	ntrajTimes,times,badTimes=ReadTimes()

	pos, vel, acc=[],[],[]
	badPos,badVel,badAcc=[],[],[]
	if readPos:
		print('ReadPos')
		pos, ntrajPos, trajlenPos, badPos=ReadObs('pos')
	if readVel:
		print('ReadVel')
		vel, ntrajVel, trajlenVel, badVel=ReadObs('vel')
	if readAcc:
		print('ReadAcc')
		acc, ntrajAcc, trajlenAcc, badAcc=ReadObs('acc')
	bad_indices=set(badPos+badVel+badAcc+badTimes)
	print('incomplete indices:',bad_indices)
	for i in bad_indices:
		if readPos and i<len(pos[i]): del pos[i] # The if i<len() is because sometimes there is nothing to remove, since bad indices come from independent lists
		if readVel and i<len(vel[i]): del vel[i]
		if readAcc and i<len(acc[i]): del acc[i]

	ntrajTimes -= len(bad_indices)
	if readPos:	
		ntrajPos -= len(bad_indices)
		checks(pos, ntrajPos, trajlenPos, 'pos')
	if readVel:	
		ntrajVel -= len(bad_indices)
		checks(vel, ntrajVel, trajlenVel, 'vel')
	if readAcc:	
		ntrajAcc -= len(bad_indices)
		checks(acc, ntrajAcc, trajlenAcc, 'acc')
		
	pos=np.array(pos,dtype=np.float64)
	acc=np.array(acc,dtype=np.float64)

	return (ntrajTimes, times, pos, vel, acc)

def skip_i(iterable, i):
    itr = iter(iterable)
    return itertools.chain(itertools.islice(itr, 0, i), itertools.islice(itr, 1, None))

def init_itimesCd(times, nt, ntCd):
	if ntCd>10:
		base=41
		fixedindices=[0,1,2,5, 9,15,22,30,base, nt-1]
		ntnew=nt-base
	else:
		return np.array([0,1]) #togliere questa riga
		fixedindices=[]
		ntnew=nt
		base=0
	ntCdnew=ntCd-len(fixedindices)
	offset=13
	fact=float((ntnew-offset)/(ntCdnew))
	otherindices=[base+offset+int(itCd*fact) for itCd in range(0,ntCdnew)]
	all_indices = np.sort(fixedindices + otherindices)
	return all_indices

class CPairPot:
	'''
	A hard-coded version of the potential used in my runs
	Potential energy between two particles.
	Potential mode is LJ xplor.
	'''
	def __init__(self, L):
		self.r_on_cutoff=np.float64(1.2)
		self.r_cutoff=np.float64(L/2-1e-10)
		self.mode='xplor'
		self.r_buff=0
		self.eps = {'AA': np.float64(1)                  , 'AB': np.float64(1.5)                , 'BB': np.float64(0.5)}
		self.sig = {'AA': np.float64(1)                  , 'AB': np.float64(0.8)                , 'BB': np.float64(0.88)}
		self.ron = {'AA': self.r_on_cutoff*self.sig['AA'], 'AB': self.r_on_cutoff*self.sig['AB'], 'BB': self.r_on_cutoff*self.sig['BB']}
		self.rcut= {'AA': self.r_cutoff   *self.sig['AA'], 'AB': self.r_cutoff   *self.sig['AB'], 'BB': self.r_cutoff   *self.sig['BB']}
		self.Vij_vectorized=np.vectorize(self.Vij, otypes=[np.float64])
		self.Vprime_vectorized=np.vectorize(self.Vprime, otypes=[np.float64])

	def __call__(self, r, pairtype):
		return self.Vij_vectorized(r, self.eps[pairtype], self.sig[pairtype], self.ron[pairtype], self.rcut[pairtype])
		# return self.Vij(r, self.eps[pairtype], self.sig[pairtype], self.ron[pairtype], self.rcut[pairtype])

	def S(self, r, r_on, r_cut):
		'''XPLOR smoothing function'''
		if r<r_on:
			return 1
		elif r>r_cut:
			return 0
		else:
			return self.Stilde(r, r_on, r_cut)
	
	def Stilde(self, r, r_on, r_cut):
		'''Stilde is the non-trivial part of S'''
		assert(r>r_on)
		assert(r<r_cut)
		rc2=r_cut*r_cut
		ro2=r_on*r_on
		r2=r*r
		term1=rc2 - r2
		term2=rc2 + 2*r2 - 3*ro2
		term3=rc2 - ro2
		return (term1*term1*term2)/(term3*term3*term3)

	
	def StildePrime(self, r, r_on, r_cut):
		'''For the derivative of S, I only consider the non-trivial part, since the rest has zero derivative'''
		assert(r>r_on)
		assert(r<r_cut)
		rc2=r_cut*r_cut
		ro2=r_on*r_on
		r2=r*r
		term1=rc2 - r2
		term2=ro2 - r2
		term3=rc2 - ro2
		return 12*r*term1*term2/(term3*term3*term3)

	def Vpair(self, r, eps, sigma, r_cut):
		'''LJ pair pure potential'''
		sigma_on_r=sigma/r
		temp=sigma_on_r*sigma_on_r
		sigma_on_r6=temp*temp*temp
		return 4*eps*(sigma_on_r6*sigma_on_r6 - sigma_on_r6)

	def VpairPrime(self, r, eps, sigma, r_cut):
		'''Derivative of the LJ pair pure potential'''
		invr=1./r
		invr2=invr*invr
		invr4=invr2*invr2
		invr6=invr4*invr2
		invr8=invr4*invr4
		s2=sigma*sigma
		s6=s2*s2*s2
		return 48*eps*s6*invr8*(0.5 - s6*invr6)

	def Vij(self, r, eps, sigma, r_on, r_cut):
		'''LJ potential with the XPLOR smoothing'''
		if r<0:
			raise ValueError('potential - S() - r=%g not admitted. Must be greater than 0.'%r)
		elif r==0: 
			return np.inf

		if r>=r_cut:
			return 0
		if(r_on<r_cut):
			return self.S(r,r_on,r_cut)*self.Vpair(r,eps, sigma,r_cut)
		elif(r_on>=r_cut):
			return self.Vpair(r,eps,sigma,r_cut)-self.Vpair(r_cut,eps,sigma,r_cut)

	def Force(self, r, pairtype):
		return -self.Vprime_vectorized(r, self.eps[pairtype], self.sig[pairtype], self.ron[pairtype], self.rcut[pairtype])


	def Vprime(self, r, eps, sigma, r_on, r_cut):
		'''Derivative of the potential with XPLOR smoothing'''
		if r<0:
			raise ValueError('potential - S() - r=%g not admitted. Must be greater than 0.'%r)
		elif r==0: 
			return -np.inf

		if r>=r_cut:
			return 0
		if self.mode=='xplor':
			if r<=r_on:
				return self.VpairPrime(r, eps, sigma, r_cut)
			else:
				return self.VpairPrime(r, eps, sigma, r_cut)*self.Stilde(r,r_on,r_cut)+(self.Vpair(r, eps, sigma, r_cut)*self.StildePrime(r, r_on, r_cut))/r
		elif self.mode=='shift' or self.mode=='no_shift':
			return self.VpairPrime(r, eps, sigma, r_cut)


def Iphi(r, normk):
	'''
	Function defined to make a sanity check in CalculateVertex().
	'''
	return integrate.quad(lambda phi: np.cos(normk*r*np.cos(phi))*np.sin(phi), 0, np.pi)[0]
Iphi_func=np.vectorize(Iphi, otypes=[np.float64])

def CalculateVertex(pos, normk):
	'''
	Calculates the vertex. The formula is written on the paper, its drafts, and in the readme file.
	'''

	t_0=time.time()

	# Number of particles of each kind
	nA = int(args.Natoms*0.8)
	nB = int(args.Natoms*0.2)
	nAA= nA*(nA-1)/2
	nBB= nB*(nB-1)/2
	nAB= nA*nB
	indicesA=np.arange(nA, dtype=int)              #First 4/5 of particles are of type A
	indicesB=np.arange(nA, args.Natoms, dtype=int) #Last  1/5 of particles are of type B
	assert(len(indicesB == nB))

	#Initialize pair potential
	U=CPairPot(L)

	#Parameters for g(r)
	dr=0.05; rMax=min(0.5*L, U.r_cutoff); nbins=int(rMax/dr)+1
	V=L*L*L; rho=args.Natoms/V

	vertex=np.ndarray(ntw)
	for itw in range(ntw):

		#Positions per type
		posA=pos[itw,0,indicesA]
		posB=pos[itw,0,indicesB]

		#Mutual distances
		distAA=med.CalculateRelativeDistances(posA, nA, L)
		distAB=np.array([med.ParticleDistPBC(posA[i], posB[j], np.array([L,L,L])) for i in range(nA) for j in range(nB)])
		distBB=med.CalculateRelativeDistances(posB, nB, L)
		distAA=np.array(list(filter(lambda x: x<=rMax, distAA)))
		distAB=np.array(list(filter(lambda x: x<=rMax, distAB)))
		distBB=np.array(list(filter(lambda x: x<=rMax, distBB)))
		nAA=len(distAA)
		nAB=len(distAB)
		nBB=len(distBB)

		rvalues,gAA=med.CalculatePairCorrelationFunction(distAA,    nA, dr=dr, rMax=rMax, number_density=nA/V)
		rvalues,gAB=med.CalculatePairCorrelationFunction(distAB, nA*nB, dr=dr, rMax=rMax, number_density=2./V)
		rvalues,gBB=med.CalculatePairCorrelationFunction(distBB,    nB, dr=dr, rMax=rMax, number_density=nB/V)

		uAA=U(rvalues,'AA')
		uAB=U(rvalues,'AB')
		uBB=U(rvalues,'BB')

		gU= (nAA*gAA[1:]*uAA[1:] + nBB*gBB[1:]*uBB[1:] + nAB*gAB[1:]*uAB[1:])/ (nAA+nBB+nAB) #The first point is nan because U(0)=infty and g(0)=0, this is why [1:]

		vertex[itw] = 4*np.pi*normk*integrate.simps( rvalues[1:]*gU*np.sin(rvalues[1:]*normk) , x=rvalues[1:])
		# vertex_check = 2*np.pi*normk*normk*integrate.simps( rvalues[1:]*rvalues[1:]*gU*Iphi_func(rvalues[1:], normk), x=rvalues[1:] )

	print('** Tempo per funzione vertex:',time.time()-t_0,'secondi')

	return vertex


def CalculateObservables(times, pos, vel, acc):
	'''
	Calculate the desired self-correlation functions
	'''
	obs={}
	ntw=max(len(pos), len(vel), len(acc)) # Need the max, because either pos,vel and acc are set to [] if when not needed
	nt=len(times)
	#permutations of the wave vector components
	print('nk = ',len(all_k))
	print('ntw = ',ntw)
	print('nt = ',nt)
	normk=np.linalg.norm(all_k[0])

	# VERTEX
	vertex=CalculateVertex(pos, normk)

	# SCATTERING FUNCTIONS

	#Array with a reduced number of times, because these observables are slow to calculate
	itimesCd=init_itimesCd(times, nt, args.ntCd)
	few_times=np.array(times)[itimesCd]
	args.ntCd=len(few_times)
	print('few_times:',few_times)
	assert(few_times[0]==0) #This is needed because of how I compressed the calculation of bkF

	#I need to loop serially over itw, otherwise too much memory is used
	t_0=time.time()
	Fkt, Fkt_all, integrand = np.ndarray((ntw, args.ntCd)), np.ndarray((ntw, args.ntCd)), np.ndarray((ntw, args.ntCd))
	for itw in range(ntw):
		print('itw:',itw)
		t_start=time.time()
		#Big matrix of displacements
		pos_few=np.array([pos[itw][it] for it in itimesCd])
		all_disps=np.array([ [ [  med.PeriodicDisplacement(pos_few[itCd][i], pos_few[0][j],L)  for j in range(args.Natoms)]for i in range(args.Natoms)] for itCd in range(args.ntCd)])
		t_alldisps=time.time()

		print('** Tempo per creare all_disps:',t_alldisps-t_start,'secondi')

		# print('np.shape(all_disps) = ',np.shape(all_disps))
		# raise KeyboardInterrupt

		Fkt_all_vec=np.cos(  np.tensordot(all_disps,all_k, axes=([3],[1]))    ).mean(axis=3)  # np.tensordot(all_disps,all_k, axes=([3],[1]))= \vec k \cdot \vec deltar
		t_medio=time.time()
		# print('** Tempo per creare Fkt_all_vec:',t_medio-t_alldisps,'secondi')

		#Calculate Self-intermediate scattering function
		Fkt[itw]=np.trace(Fkt_all_vec, axis1=1, axis2=2)/(args.Natoms)
		#Calculate Collective scattering function
		Fkt_all[itw]=Fkt_all_vec.sum(axis=(1,2))

		# A single entry for the integrand in dk
		integrand[itw]=normk*normk*vertex[itw]*vertex[itw]*(Fkt[itw]*Fkt_all[itw])
		t_integrand=time.time()
		# print('** Time to calculate the general scattering function:', int((t_integrand-t_medio)),' seconds')

	print('* Time to calculate the general scattering function:', int((t_integrand-t_0)/60),'minutes')


	# Save correlation functions on an array
	obs['few_times'] = few_times
	obs['Fkt'      ] = {'mean': np.mean(      Fkt, axis=0), 'err': sem(      Fkt, axis=0), 'blocks':       Fkt}
	obs['Fkt_all'  ] = {'mean': np.mean(  Fkt_all, axis=0), 'err': sem(  Fkt_all, axis=0), 'blocks':   Fkt_all}
	obs['vertex'   ] = {'mean': np.mean(   vertex, axis=0), 'err': sem(   vertex, axis=0), 'blocks':    vertex}
	obs['integrand'] = {'mean': np.mean(integrand, axis=0), 'err': sem(integrand, axis=0), 'blocks': integrand}

	return normk, obs


def WriteCorrelations(k):
	'''
	Writes to file the measured correlations, both in binary (for further analysis) and in text (for plots)
	'''
	print('WriteCorrelations')
	if args.thermostat=='*': 
		args.thermostat=''
	np.savetxt('vertexObs_{}_k{}.txt'.format(args.thermostat,k), 
		np.column_stack((
			obs['few_times'], 
			args.n[0]*np.ones(args.ntCd),args.n[1]*np.ones(args.ntCd),args.n[2]*np.ones(args.ntCd),
			k*np.ones(args.ntCd),
			obs['Fkt']['mean'], obs['Fkt']['err'],
			obs['Fkt_all']['mean'], obs['Fkt_all']['err'],
			obs['vertex']['mean']*np.ones(args.ntCd), obs['vertex']['err']*np.ones(args.ntCd),
			obs['integrand']['mean'], obs['integrand']['err']
			)), 
		fmt=['%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g'], 
		header='time n1 n2 n3 k Fkt(self) errFktSelf Fkt(collective) errFktColl vertex errVertex integrand errIntegrand')
	np.save('vertexObs_{}_k{}.npy'.format(args.thermostat,k), obs)
	return




if __name__=='__main__':
	dawn=time.time()
	(ntw,times,pos,vel,acc)=ReadAll()
	normk,obs=CalculateObservables(times,pos,vel,acc)
	WriteCorrelations(normk)
	dusk=time.time()
	print('*** Total simulation time:',(dusk-dawn)/60,'minutes')



