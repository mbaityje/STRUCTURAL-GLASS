#!/usr/bin/env python

import sys, argparse, time, itertools, re, numpy as np
from matplotlib import pyplot as plt
from scipy.stats import sem
from scipy import integrate, interpolate
from glob import glob
from math import isclose

def k_from_name(name):
	'''
	Extract the k from the filename, and double check the value on the corresponding .txt file
	'''
	remove_init=re.sub(r'^(.*)k',"",name)
	k=np.float64(re.sub(r'.npy$',"",remove_init))
	
	#Double check with the txt file
	nametxt=re.sub(r'.npy','.txt',name)
	ktxt=np.loadtxt(nametxt)[0,1:5] #ktxt=(nx, ny, nz, |k|)
	#Check that the modulus is correct
	if np.any(np.isclose(ktxt[3],k)) == False:
		raise ValueError('The k read in the txt file is different from the one on the filename')
	#Check that the modulus and the quantum numbers agree
	kvec=ktxt[:3]*2*np.pi/args.L
	if np.any(np.isclose(np.linalg.norm(kvec),k)) == False:
		print('The k read in line 4 of the txt file is not the modulus of the previous three lines')
		print('kvec = {}, k = {}'.format(kvec,k))
		raise ValueError('The k read in line 4 of the txt file is not the modulus of the previous three lines')
	return k, kvec

def kvec_from_nametxt(nametxt):
	return np.loadtxt(all_namestxt[0])[0][1:4]

parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--L', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('--thermostat', required=False, default='NVT', help='thermostat')
parser.add_argument('--dt', required=False, default=0.0025, help='thermostat')
parser.add_argument('--kmax', type=np.float, required=False, default=None, help='Maximum value of k over which to integrate')
parser.add_argument('--showplots', action='store_true', help='If flag is activated, makes some plots. Each plot stops the program flow, which must be resumed manually.')
parser.add_argument('--fuchs', action='store_true', help='If flag is activated, uses the Fuchs approach, from Fuchs&Gotze(1999) and Balucani&Zoppi (pg.186).')
args = parser.parse_args()

interpkind='quadratic'

##############################
# Read observables for all k #
##############################

names='./vertexObs_{}_k*.npy'.format(args.thermostat)
all_names=sorted(glob(names)) #sorted gives alphabetical sorting because it's strings, so must add a couple of lines
klist=sorted([k_from_name(name)[0] for name in all_names], key=np.float64)
kmin=klist[0]
nk=len(klist)

all_names   =['./vertexObs_{}_k{}.npy'.format(args.thermostat,str(k)) for k in klist]
all_namestxt=['./vertexObs_{}_k{}.txt'.format(args.thermostat,str(k)) for k in klist]
obs = [ (k_from_name(name)[0], k_from_name(name)[1], np.load(name).tolist()) for name in all_names] # contains all the observables for each k

#Check that we stored the right k: compare components and modulus
for ik in range(nk):
	if not isclose(obs[ik][0], np.linalg.norm(obs[ik][1])):
		raise ValueError('Conflict between |k| and kvec')

nt=len(obs[0][2]['few_times'])
rho=args.Natoms/(args.L*args.L*args.L)
times=obs[0][2]['few_times']*args.dt


assert (nk==len(obs) ) #Make sure all files were read
assert( np.all( [len(obs[ik][2]['few_times'])==nt for ik in range(nk)])==True ) #Make sure all datafiles have the same number of times

###########################
# Static structure factor #
###########################

Sk=np.array([ obs[ik][2]['Fkt_all']['mean'][0] for ik in range(nk)])/args.Natoms
Sk_err=np.array([ obs[ik][2]['Fkt_all']['err'][0] for ik in range(nk)])/args.Natoms
Sk_interp=interpolate.interp1d(klist, Sk, kind=interpkind)
if args.showplots:
	#Plot S(k)
	plt.title('Static Structure Factor')
	plt.ylabel('$S(k)$')
	plt.xlabel('$k$')
	plt.errorbar(klist, Sk, yerr=Sk_err, color='red')
	ascisse=np.arange(klist[0],klist[-1],0.1)
	plt.plot(ascisse, Sk_interp(ascisse), color='red')
	plt.plot(klist, np.ones(nk), linestyle='dashed', color='black', linewidth=0.1)
	plt.show()


def ComputeK(obs, it, kmax=None, fuchs=False):
	'''
	Calculate the MCT-kind memory function. 
	Allows for two possible method. One is the MCT that we detail in the appendix of one of the versions,
	and the other one is from a paper by Fuchs&Gotze.
	The flag `fuchs` switches from between the two.
	'''
	klist = np.array([obs[ik][0] for ik in range(nk)])
	kvec  = np.array([obs[ik][1] for ik in range(nk)])
	FF    = np.array( [ obs[ik][2]['Fkt']['mean'][it] for ik in range(nk)])*np.array( [ obs[ik][2]['Fkt_all']['mean'][it] for ik in range(nk)]) #Fkt*Fkt_all
	ck = (1-1./Sk)/rho
	ck2=ck*ck

	if fuchs==False:
		integrand     = np.array([ obs[ik][2]['integrand']['mean'][it] for ik in range(nk)])

		#In case we want to recalculate it to make sure that the integrand saved on disk is correct:
		# integrand     = np.array( [ obs[ik][2]['Fkt']    ['mean'][it] for ik in range(nk)])*\
		# 				np.array( [ obs[ik][2]['Fkt_all']['mean'][it] for ik in range(nk)])*\
		# 				np.square([ obs[ik][2]['vertex' ]['mean']     for ik in range(nk)])*\
		# 				klist*klist

		integrand_interp=interpolate.interp1d(klist, integrand, kind=interpkind)
		if kmax==None: kmax=klist[-1]
		output=integrate.quad(integrand_interp, klist[0], kmax, maxp1=100)[0] / (2*np.pi*np.pi*rho*args.temperature)
	elif fuchs==True:
		# The vertex in this case is different, vertex=ck*k_alpha:
		prefactor=args.temperature*rho/(2*np.pi*np.pi)
		output=0
		for idim in range(3):
			integrand=klist*klist*kvec[:,idim]*kvec[:,idim]*ck2*FF
			integrand_interp=interpolate.interp1d(klist, integrand, kind=interpkind)
			output += prefactor*integrate.quad(integrand_interp, klist[0], kmax, maxp1=100)[0]
		output/=3

	return output



unbounded=False
if args.fuchs==False and args.kmax==None: # Compare with correlators at time zero, to use the best kmax. For fuchs=True we don't have a cutoff.
	F0=np.load('../CFF_NVT.npy').item()['mean'][0]
	K0=F0/args.temperature
	rtol=0.02
	err=np.inf
	kstep=0.01

	for k in np.arange(kmin, klist[-1], kstep):
		Ktemp=ComputeK(obs, 0, kmax=k, fuchs=args.fuchs)
		if np.abs(K0-Ktemp)>err:
			raise StopIteration('Unable to get close to K0. Try reducing the step size kstep')
		err=np.abs(K0-Ktemp)
		print('k:%.2f\tK0: %g Ktemp=%g err=%g'%(k,K0,Ktemp,err))
		if err < rtol*K0:
			args.kmax=k
			break
elif args.kmax==None or args.kmax>klist[-1]:
	args.kmax=klist[-1]
	unbounded=True

K=np.zeros(nt)
for it in range(nt):
	K[it] = ComputeK(obs, it, args.kmax, fuchs=args.fuchs)

#Save Memory function
# name = '../Kvertex_{}_kmax{}.txt'.format(args.thermostat,args.kmax) if (not args.fuchs) else '../KvertexFuchs_{}_kmax{}.txt'.format(args.thermostat,args.kmax)


name= '../Kvertex'+ ('Fuchs_' if args.fuchs else '_') + args.thermostat +   ('.txt' if unbounded else '_kmax%g.txt'%args.kmax)

np.savetxt(name,
	np.column_stack((times, K, K/K[0])),
	fmt=['%g','%.14g','%.14g'],
	header='time Kvertex Kvertex/Kvertex[0]'
	)

#Save static structure factor
np.savetxt('../Sk_{}.txt'.format(args.thermostat),
	np.column_stack((klist, Sk, Sk_err)),
	fmt=['%g','%.14g','%.14g'],
	header='time S(k) errS(k)'
	)


if args.showplots:

	#Plot memory function
	plt.title('MEMORY FUNCTION')
	plt.plot(times,K/K[0])
	plt.xscale('log')
	plt.show()

