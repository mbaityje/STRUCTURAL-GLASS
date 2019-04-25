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
	np.loadtxt(nametxt)[:,4]
	if np.any(np.isclose(np.loadtxt(nametxt)[:,4],k)) == False:
		raise ValueError('The k read in the txt file is different from the one on the filename')
	return k

parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--L', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('--thermostat', required=False, default='NVT', help='thermostat')
parser.add_argument('--dt', required=False, default=0.0025, help='thermostat')
parser.add_argument('--kmax', type=np.float, required=False, default=None, help='Maximum value of k over which to integrate')
parser.add_argument('--showplots', action='store_true', help='If flag is activated, makes some plots. Each plot stops the program flow, which must be resumed manually.')
args = parser.parse_args()

##############################
# Read observables for all k #
##############################

names='./vertexObs_{}_k*.npy'.format(args.thermostat)
all_names=sorted(glob(names)) #sorted gives alphabetical sorting because it's strings, so must add a couple of lines
klist=sorted([k_from_name(name) for name in all_names], key=np.float64)
kmin=klist[0]
nk=len(klist)

all_names=['./vertexObs_{}_k{}.npy'.format(args.thermostat,str(k)) for k in klist]
obs = [ (k_from_name(name), np.load(name).tolist()) for name in all_names] # contains all the observables for each k
nt=len(obs[0][1]['few_times'])
rho=args.Natoms/(args.L*args.L*args.L)
times=obs[0][1]['few_times']*args.dt

assert (nk==len(obs) ) #Make sure all files were read
assert( np.all( [len(obs[ik][1]['few_times'])==nt for ik in range(nk)])==True ) #Make sure all datafiles have the same number of times

###########################
# Static structure factor #
###########################

Sk=np.array([ obs[ik][1]['Fkt_all']['mean'][0] for ik in range(nk)])/args.Natoms
Sk_err=np.array([ obs[ik][1]['Fkt_all']['err'][0] for ik in range(nk)])/args.Natoms
Sk_interp=interpolate.interp1d(klist, Sk, kind='cubic')
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


def ComputeK(klist, obs, it, kmax=None):
	integrand     = np.array([ obs[ik][1]['integrand']['mean'][it] for ik in range(nk)])

	#In case we want to recalculate it to make sure that the integrand saved on disk is correct:
	# integrand     = np.array( [ obs[ik][1]['Fkt']    ['mean'][it] for ik in range(nk)])*\
	# 				np.array( [ obs[ik][1]['Fkt_all']['mean'][it] for ik in range(nk)])*\
	# 				np.square([ obs[ik][1]['vertex' ]['mean']     for ik in range(nk)])*\
	# 				klist*klist

	integrand_interp=interpolate.interp1d(klist, integrand, kind='quadratic')
	if kmax==None: kmax=klist[-1]
	return integrate.quad(integrand_interp, klist[0], kmax)[0] / (2*np.pi*np.pi*rho*args.temperature)



if args.kmax==None: # Compare with correlators at time zero, to use the best kmax
	F0=np.load('../CFF_NVT.npy').item()['mean'][0]
	K0=F0/args.temperature
	rtol=0.02
	err=np.inf
	kstep=0.01

	for k in np.arange(kmin, klist[-1], kstep):
		Ktemp=ComputeK(klist, obs, 0, kmax=k)
		if np.abs(K0-Ktemp)>err:
			raise StopIteration('Unable to get close to K0. Reduce step size kstep')
		err=np.abs(K0-Ktemp)
		print('k:%.2f\tK0: %g Ktemp=%g err=%g'%(k,K0,Ktemp,err))
		if err < rtol*K0:
			args.kmax=k
			break
elif args.kmax>klist[-1]:
	print('WARNING: args.kmax=%g is larger than the largest simulated k, ktop=%g, so we set kmax=ktop'%(args.kmax,klist[-1]))
	args.kmax=klist[-1]

K=np.zeros(nt)
for it in range(nt):
	K[it] = ComputeK(klist, obs, it, args.kmax)
	

#Save Memory function
np.savetxt('../Kvertex_{}.txt'.format(args.thermostat),
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

