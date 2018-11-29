#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
import lib.module_measurements as med
from matplotlib import pyplot as plt
from lib.beylkin import Beylkin


def FitFunction(x, y, ncoef=10):
	'''Fit function as a sum of exponential functions'''
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function


# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--box_size', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='time step')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--maxtime', type=float, required=False, default='-1.0', help='truncate the time axis at maxtime')
parser.add_argument('--tstar', type=float, required=False, default=None, help='no integration for t>tstar')
parser.add_argument('--shiftCFP', action='store_true', help='if invoked, imposes CFP[0]=0')
parser.add_argument('--softening', action='store_true', help='if invoked, soften kernels in order to have less signal where signal is crap')
parser.add_argument('--lin', action='store_true', help='calculate noise correlation on linear grid')
parser.add_argument('--linsc', action='store_true', help='calculate noise correlation self-consistently on linear grid')
parser.add_argument('--normalsc', action='store_true', help='calculate noise correlation self-consistently (on generic grid)')
parser.add_argument('--fits', action='store_true', help='do fits instead of interpolations')

args = parser.parse_args()

print(sys.argv)
print('N:',args.Natoms)
print('L:',args.box_size)
print('T:',args.temperature)
print('dt:',args.dt)
print('thermostat:',args.thermostat)
L=args.box_size

if not args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
if not args.box_size>0: raise ValueError('L ({}) must be positive'.format(args.box_size))
if not args.dt>0: raise ValueError('dt ({}) must be positive'.format(args.dt))
if not args.dt<0.1: raise ValueError('I do not believe dt is that big ({})'.format(args.dt))
if (not args.tstar==None) and args.tstar<=0: raise ValueError('tstar should be positive ({})'.format(args.tstar))
if args.temperature<0:
	raise ValueError('Temperature='+str(args.temperature)+'. Must be positive. Aborting.')
else:
	invT=1./args.temperature



#READ KERNELS
CFF=np.load('CFF_'+str(args.thermostat)+'.npy') # access as CFP.item()['mean']
CFP=np.load('CFP_'+str(args.thermostat)+'.npy')
times=np.load('times_'+str(args.thermostat)+'.npy')
nt=len(times)
if not (len(CFP.item()['mean'])==nt and len(CFF.item()['mean'])==nt):
	raise IndexError('The number of correlations ('+str(len(CFP.item()['mean']))+' or '+str(len(CFF.item()['mean']))+') is inconsistent with the number of times('+str(nt)+')')

invT=np.float64(1./args.temperature)
dtinvT=args.dt*invT

# PREPARE KERNELS
#Consider truncation of the time axis
if args.maxtime>0:
	times=np.array([t for t in times if t<args.maxtime])
	nt=len(times)
	CFF.item()['mean']=CFF.item()['mean'][:nt]
	CFF.item()['err' ]=CFF.item()['err' ][:nt]
	CFP.item()['mean']=CFP.item()['mean'][:nt]
	CFP.item()['err' ]=CFP.item()['err' ][:nt]
print('Integration interval: [{},{}]'.format(times[0],times[nt-1]))

if args.shiftCFP:
	CFP.item()['mean'][0]-=CFP.item()['mean'][0] #Do not shift the whole curve, or there will be a bias in the integral!!!

#Softening
if args.softening:
	raise NotImplementedError('No softening yet, gotta think a rigorous method')
	isoft=int(nt*.2)
	print('isoft = ',isoft)
	CFP.item()['mean'][isoft:]*=[np.exp(-i/100) for i in range(0,nt-isoft)]
	CFF.item()['mean'][isoft:]*=[np.exp(-i/100) for i in range(0,nt-isoft)]


def NoiseCorrLinear(times, CFF, CFP, dt=0.0025):
	'''
	Calculate noise correlation function on a linear grid, with the non-selfconsistent method.
	Input is data on any grid, and the linear one is generated with time step dt.
	'''
	#Create CFF and CFP on a linear grid
	def trapezeLinear(i,istar=None):
		temp=0.5*kernel[i]*f[0]
		i_ini=1 if istar==None else max(1,i-istar)
		for j in range(i_ini, i-1):
			temp+=kernel[i-j]*f[j]
		return g[i]+dtinvT*temp

	#Interpolate CFF and CFP
	interpCFF=FitFunction(times, CFF, ncoef=10) if args.fits else interp1d(times, CFF, kind='cubic')
	interpCFP=FitFunction(times, CFP, ncoef=10) if args.fits else interp1d(times, CFP, kind='cubic')


	#Create linear grid and observables on it
	lineargrid=np.arange(0, times[-1], dt)
	ntlinear=len(lineargrid)
	f=np.zeros(ntlinear)
	g=interpCFF(lineargrid)
	kernel=interpCFP(lineargrid)

	if args.tstar==None:
		istar=None #for i>istar the kernels have no signal
	else:
		for i in range(len(lineargrid)):
			if lineargrid[i]>args.tstar:
				istar=i
				break

	print('istar = ', istar)
	if not istar==None: print('t[istar] = ',lineargrid[istar])

	#Calculate the Volterra solution
	f[0]=g[0]
	for i in range(1, ntlinear):
		print('\rtrapezeLinear iteration',i, end='')
		f[i]=trapezeLinear(i,istar)
	print('')
	return lineargrid, kernel, g, f


from scipy.integrate import simps,romb
def NoiseCorrSelfConsistentLinear(lineargrid, linearCFF, linearCFP, Kold=None, scheme='trapeze', maxiter=1000):
	'''
	Calculate the noise correlation function through the self-consistent equation on a linear grid.
	Input must be functions on a linear scale.
	Kold: initial guess for the correlation function
	'''
	print('Checking with self-consistent formulation on linear grid')

	def integral(myCFP, myKold, myn, dt, scheme='rectangles'):
		'''Integral for the self-consistent calculation'''
		if n==0: return 0
		integrand=np.array([myCFP[myn-i]*myKold[i] for i in range(myn)])

		if scheme   == 'rectangles':
			output=integrand.sum()*dt
		elif scheme == 'trapeze':
			output=np.trapz(integrand,dx=dt)
		elif scheme == 'simpson':
			output=simps(integrand,dx=dt, even='first')
		else:
			raise NotImplementedError('Integration scheme must be rectangles, trapeze or simpson. '+scheme+' is not a valid option.')
		return output

	nt=len(lineargrid)
	dt=lineargrid[1]-lineargrid[0]
	assert(np.isclose(dt,lineargrid[-1]-lineargrid[-2])) # Small check that grid is linear
	if Kold is None:
		Kold=linearCFF
	Knew=np.zeros(nt)

	print('len(Kold): ',len(Kold))
	print('nt: ',nt)
	assert(nt==len(linearCFF)) # Check that grid and data are consistent
	assert(nt==len(linearCFP))
	assert(nt==len(Kold))

	for it in range(maxiter):
		for n in range(nt):
			Knew[n] = linearCFF[n] + invT * integral(linearCFP, Kold, n, args.dt, scheme=scheme)
		err=np.max(np.abs(Knew-Kold))
		print("it:",it," err = ",err)
		if err<1e-8: break
		Kold[:]=Knew[:]
	return Knew


def NoiseCorr(times, CFF, CFP):
	'''
	Calculate noise correlation function on a GENERIC grid, with the non-selfconsistent method.
	'''
	nt=len(times)
	assert(nt==len(CFF) and nt==len(CFP))
	f=np.zeros(nt)
	g=CFF
	def kernel(i,j):
		if i==j: return 0
		return interpCFP(times[i]-times[j])

	#Interpolate CFF and CFP
	# interpCFF= FitFunction(times, CFF, ncoef=10) if args.fits else  interp1d(times, CFF, kind='cubic')
	interpCFP= FitFunction(times, CFP, ncoef=20) if args.fits else  interp1d(times, CFP, kind='cubic')


	#eventualmente mettere istar

	#Create CFF and CFP on a generic grid
	def trapeze(i):
		temp=kernel(i,i-1)*f[i-1]*(times[i]-times[i-1])
		for j in range(1, i-1):
			temp+=(kernel(i,j)*f[j]+kernel(i,j-1)*f[j-1])*(times[j]-times[j-1])
		return g[i]+0.5*invT*temp

	#Calculate the Volterra solution
	f[0]=g[0]
	for i in range(1, nt):
		print('\rtrapeze iteration',i, end='')
		f[i]=trapeze(i)
		if f[i]>f[0]:
			raise ValueError('NoiseCorr is giving unphysical values (it is growing)')
	print('')
	return f



from scipy.integrate import quad
def NoiseCorrSelfConsistent(times, CFF, CFP, Kold=None, maxiter=1000, tstar=None):
	'''
	Calculate the noise correlation function through the self-consistent equation.
	Kold: initial guess for the correlation function
	'''

	print('Checking with self-consistent formulation on generic grid')

	nt=len(times)
	if Kold is None:
		Kold=CFF
	Knew=np.zeros(nt)
	assert(nt==len(CFF)) # Check that grid and data are consistent
	assert(nt==len(CFP))
	assert(nt==len(Kold))
	if not tstar==None: 
		istar=np.where(times>args.tstar)[0][0]
		maxn=istar # redundant, but notation is important too :P
	else: maxn=nt

	interpCFP  = FitFunction(times, CFP , ncoef=10) if args.fits else interp1d(times, CFP , kind='cubic')
	interpKold = FitFunction(times, Kold, ncoef=10) if args.fits else interp1d(times, Kold, kind='cubic')

	print('maxn = ',maxn)
	print('nt = ',nt)

	for it in range(maxiter):
		for n in range(maxn):
			print('\r\tn:',n,end='')
			t=times[n]
			temp=quad(lambda u: interpCFP(t-u)*interpKold(u), 0, t, limit=20, maxp1=20, limlst=20)[0]
			Knew[n] = CFF[n] + invT * temp
		print('')
		Kold=interpKold(times)
		err=np.max(np.abs(Knew-Kold))
		print("it:",it," err = ",err)
		if err<1e-8: 
			break
		interpKold=FitFunction(times, Knew, ncoef=10) if args.fits else interp1d(times, Knew, kind='cubic')
	return Knew







#                         #
# HERE STARTS THE PROGRAM #
#                         #

# Calculate noise correlation function on the measurement grid
f = NoiseCorr(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'])
np.savetxt('noisecorr_{}.txt'.format(args.thermostat), 
			np.column_stack((times, f, f/f[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time K K/K[0]')

if args.normalsc:
	fMakesSense = True if all([item<=f[0] for item in f]) else False
	print('fMakesSense:', fMakesSense)
	fsc=NoiseCorrSelfConsistent(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'], Kold=f if fMakesSense else None, tstar=args.tstar)
	np.savetxt('noisecorr-selfcon_{}.txt'.format(args.thermostat), 
			np.column_stack((times, fsc, fsc/fsc[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time K K/K[0]')

# Calculate noise correlation function on linear grid
if args.lin:
	(lineartimes,linearCFP, linearCFF, flin) = NoiseCorrLinear(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'], dt=args.dt)
	np.savetxt('noisecorrlinear_{}.txt'.format(args.thermostat), 
				np.column_stack((lineartimes, flin, flin/flin[0])), 
				fmt=['%.14g','%.14g','%.14g'], 
				header='time K K/K[0]')

# Calculate noise correlation function selfconsistently on linear grid
if args.linsc:
	flinMakesSense = True if all([item<=flin[0] for item in flin]) else False
	flinsc = NoiseCorrSelfConsistentLinear(lineartimes, linearCFF, linearCFP, Kold=flin if flinMakesSense else None)
	np.savetxt('noisecorrlinear-selfcon_{}.txt'.format(args.thermostat), 
				np.column_stack((lineartimes, flinsc, flinsc/flinsc[0])), 
				fmt=['%.14g','%.14g','%.14g'], 
				header='time K K/K[0]')













