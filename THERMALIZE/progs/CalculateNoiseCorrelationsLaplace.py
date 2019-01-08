#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
from scipy.integrate import simps,romb,quad
import lib.module_measurements as med
from matplotlib import pyplot as plt
from lib.beylkin import Beylkin
from lib.module_timelists import ListaLogaritmica

def NoiseCorr(times, CFF, CFP):
	'''
	Calculate noise correlation function on a GENERIC grid, with the non-selfconsistent method.
	'''
	nt = len(times)
	assert(nt==len(CFF) and nt==len(CFP))
	f=np.zeros(nt)
	g=CFF
	def kernel(i,j):
		if i==j: return 0
		return interpCFP(times[i]-times[j])

	#Interpolate CFF and CFP
	# interpCFF= FitFunction(times, CFF, ncoef=10) if args.fits else  interp1d(times, CFF, kind='cubic')
	interpCFP= FitFunction(times, CFP, ncoef=20) if args.fits else  interp1d(times, CFP, kind='cubic', assume_sorted=True)

	#eventualmente mettere istar

	#Create CFF and CFP on a generic grid
	def trapeze(i):
		temp=kernel(i,i-1)*f[i-1]*(times[i]-times[i-1])
		for j in range(1, i-1):
			temp+=(kernel(i,j)*f[j]+kernel(i,j-1)*f[j-1])*(times[j]-times[j-1])
		# print('temp=',temp)
		return g[i]+0.5*invT*temp

	#Calculate the Volterra solution
	f[0]=g[0]
	for i in range(1, nt):
		print('\rtrapeze iteration',i, end='')
		f[i]=trapeze(i)
		# if f[i]>f[0]:
			# raise ValueError('NoiseCorr is giving unphysical values (it is growing)')
	print('')
	return f

def NoiseCorrSelfConsistent(times, CFF, CFP, Kold=None, maxiter=1000, tstar=None):
	'''
	Calculate the noise correlation function through the self-consistent equation.
	Input must be functions on a linear scale.
	Kold: initial guess for the correlation function
	'''

	print('Kold:',Kold)

	nt=len(times)
	if Kold is None:
		Kold=CFF
	Knew=np.zeros(nt)
	assert(nt==len(CFF)) # Check that grid and data are consistent
	assert(nt==len(CFP))
	assert(nt==len(Kold))
	if not tstar==None: 
		istar=np.where(times>tstar)[0][0]
		maxn=istar # redundant, but notation is important too :P
	else: maxn=nt

	interpCFP  = FitFunction(times, CFP , ncoef=10) if args.fits else interp1d(times, CFP , kind='cubic', assume_sorted=True)
	interpKold = FitFunction(times, Kold, ncoef=10) if args.fits else interp1d(times, Kold, kind='cubic', assume_sorted=True)

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
		plt.semilogx(times, Knew)
		plt.show()
		if err<1e-8: 
			break
		interpKold=FitFunction(times, Knew, ncoef=10) if args.fits else interp1d(times, Knew, kind='cubic', assume_sorted=True)
	return Knew

# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
parser.add_argument('-L','--box_size', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='time step')
parser.add_argument('-M','--M', type=int, required=False, default=3, help='Number of particles')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--maxtime', type=float, required=False, default='-1.0', help='truncate the time axis at maxtime')
parser.add_argument('--tstar', type=float, required=False, default=None, help='no integration for t>tstar')
parser.add_argument('--shiftCFP', action='store_true', help='if invoked, imposes CFP[0]=0')
parser.add_argument('--softening', action='store_true', help='if invoked, soften kernels in order to have less signal where signal is crap')
parser.add_argument('--fits', action='store_true', help='do fits instead of interpolations')
parser.add_argument('--kind', required=False, choices=['interp','interp_lin','fit','combined'], default='interp', help='thermostat')
parser.add_argument('--ncoef', type=int, required=False, default=10, help='number of coefficients for fits')
parser.add_argument('--showplots', action='store_true', help='if invoked, shows plots during the calculation (this stops the calculation until the user closes the figure)')


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
CPP=np.load('CPP_'+str(args.thermostat)+'.npy')
times=np.load('times_'+str(args.thermostat)+'.npy')

# I want to know the position of t=1 LJ unit
it_one = np.where(times>=1)[0][0]

nt=len(times)
if not (len(CFP.item()['mean'])==nt and len(CFF.item()['mean'])==nt and len(CPP.item()['mean'])==nt):
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
	CPP.item()['mean']=CPP.item()['mean'][:nt]
	CPP.item()['err' ]=CPP.item()['err' ][:nt]
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
	CPP.item()['mean'][isoft:]*=[np.exp(-i/100) for i in range(0,nt-isoft)]


def FitFunction(x, y, ncoef=10):
	'''Fit function as a sum of exponential functions'''
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function

def LaplaceSlow(f, p, tmax=np.inf):
	return quad(lambda t: f(t)*np.exp(-t*p), 0, tmax)[0]

def GaverStehfest(ftilde, t, M):
	if t==0:
		print('time t=%g is too small'%t)
		return np.nan
	fact=np.math.factorial
	l2ont=np.log(2)/t
	def omega(k):
		sign = 1 if (M+k)%2==0 else -1
		top=min(k,M)
		bottom=int((k+1)/2)
		val=0
		for j in range(bottom, top+1):
			val += j**(M) * fact(2*j) / ( fact(M-j) * fact(j) * fact(j-1) *fact(k-j) * fact(2*j-k))
		return sign*val
	out=0
	for k in range(1, 2*M+1):
		# print('k*log(2)/t = ',k*l2ont)
		out += omega(k) * ftilde(k*l2ont)
	return out*l2ont



def TransformAntitransform(x, y, M, ncoef=10, showplots=False, kind='interp'):
	'''
	Takes an array y(x), Laplace transforms it and antitransforms it.
	This tells us if we are in the right parameter range.
	'''

	#The array of laplace space impulses
	pmin=np.log(2)/x[-1]-1e-8
	pmax=2*M*np.log(2)/(x[1]-x[0])
	deltap=2*M*np.log(2)/(x[-1]-x[-2])
	print('pmin:',pmin,'pmax:',pmax, 'deltap:',deltap)
	p=np.arange(pmin, pmax+deltap, deltap)
	p=ListaLogaritmica(pmin, pmax+1e-3, 500)

	# func=FuncFromArray(x, y, kind=kind, ncoef=ncoef)
	func=FuncFromArray(x, y, kind=kind, ncoef=ncoef)


	if showplots:
		plt.title('Pure function')
		plt.semilogx(x,y, label='data')
		plt.semilogx(x,func(x), label=kind+' interpolation')
		plt.legend()
		plt.show()

	print('Transforming...')
	# transform=np.array([LaplaceSlow(func, xp, tmax=1) for xp in x])
	transform=np.ndarray(len(p))
	for ip in range(len(p)):
		pp=p[ip]
		tmax=x[-1] if (kind=='interp' or kind=='interp_lin') else np.inf
		val=LaplaceSlow(func, pp, tmax=tmax)
		print('\rpp={},  Lf={}'.format(pp, val), end='')
		transform[ip]=val
	print('')

	func_transform=FuncFromArray(p, transform, kind='interp', ncoef=ncoef)

	if showplots:
		plt.title('Transform')
		plt.semilogx(p, transform, label='data')
		plt.semilogx(p, func_transform(p), label=kind+' interpolation')
		plt.legend()
		plt.show()

	print('Anti-Transforming...')
	# antitransform=np.array([GaverStehfest(func_transform, xp) for xp in x])
	antitransform=np.ndarray(len(x))
	for ix in range(len(x)):
		xp=x[ix]
		val=GaverStehfest(func_transform, xp, M=args.M)
		print('\rxp=%g,  LinvLf=%g'%(xp,val), end='')
		antitransform[ix]=val
	print('')

	if showplots:
		plt.title('Pure function')
		plt.semilogx(x,y, label='data')
		plt.semilogx(x,func(x), label='fit')
		plt.semilogx(x,antitransform, label='antitransform')
		plt.legend()
		plt.show()

	return p,transform, antitransform

def FuncFromArray(x, y, kind='interp', ncoef=10):
	if kind == 'interp':
		return interp1d(x, y, kind='cubic', assume_sorted=True)
	elif kind == 'interp_lin':
		return interp1d(x, y, kind='slinear', assume_sorted=True)
	elif kind == 'fit':
		return FitFunction(x,y, ncoef=ncoef)
	elif kind == 'combined':
		tempfunc=interp1d(x, y, kind='cubic', assume_sorted=True)
		mydelta=(x[it_one]-x[0])/1000
		xlinear=np.arange(x[0], x[it_one]+mydelta, mydelta)
		return FitFunction(xlinear,tempfunc(xlinear), ncoef=ncoef)
	else:
		raise NotImplementedError('FuncFromArray: not implemented option')




# ELIMINARE (prova per vedere se i fit vengono bene)

# sys.exit('Fine Prova')
# # b = Beylkin(decaying=True)
# datapoints=CPP.item()['mean']
# # b.driver_load(times, datafunc, 10)
# datafunc=FuncFromArray(times, datapoints, kind='fit', ncoef=10)
# plt.semilogx(times, datapoints, label='data')
# t=np.arange(0,times[-1],0.001)
# plt.semilogx(times,datafunc(times), label='fit')
# plt.legend()
# plt.show()








#                         #
# HERE STARTS THE PROGRAM #
#                         #

cpp=np.copy(CPP.item()['mean'])
tlist=np.copy(times)

# pcff,Lcff,LLcff=TransformAntitransform(times, CFF.item()['mean'], ncoef=10, showplots=True, kind=args.kind)
# pcfp,Lcfp,LLcfp=TransformAntitransform(times, CFP.item()['mean'], ncoef=10, showplots=True, kind=args.kind)
pcpp,Lcpp,LLcpp=TransformAntitransform(times, cpp, M=args.M, ncoef=args.ncoef, showplots=args.showplots, kind=args.kind)




LK = np.ndarray(len(Lcpp))
for ip in range(len(Lcpp)):
	LK[ip] = ( args.temperature - pcpp[ip]*Lcpp[ip] ) / Lcpp[ip]

LLK=np.ndarray(len(tlist))
func_LK=FuncFromArray(pcpp, LK, kind='interp', ncoef=args.ncoef)

for it in range(len(tlist)):
	tp=tlist[it]
	val=GaverStehfest(func_LK, tp, M=args.M)
	print('\rxp=%g,  LinvLf=%g'%(tp,val), end='')
	LLK[it]=val
print('')


Kvolterra = NoiseCorr(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'])/args.temperature

if args.showplots:
        plt.title('Noise Correlation')
        plt.xlabel('$t$')
        plt.ylabel('$\mathcal{K}(t)$')
        plt.semilogx(times, LLK, label='$\\frac{kT-sC^P(s)}{C^P(s)}$')
        plt.semilogx(times, Kvolterra, label='K(t) [Volterra]')
        plt.legend()
        plt.show()

# Combine Kvolterra with LLK
tmin=0.05
if args.temperature==5.0:
	tmin=0.027
if args.temperature==2.0:
	tmin=0.042
elif args.temperature==1.0:
	tmin=0.043
elif args.temperature==0.8:
	tmin=0.04
elif args.temperature==0.7:
	tmin=0.05
elif args.temperature==0.6:
	tmin=0.05
elif args.temperature==0.55:
	tmin=0.05
itmin=np.where(times>=tmin)[0][0]

Kcombine=np.copy(Kvolterra)
Kcombine[itmin:]=LLK[itmin:]


if args.showplots:
        plt.title('Noise Correlation')
        plt.semilogx(times, Kcombine, label='K(t) [Combine]')
        plt.legend()
        plt.show()

#Same Kcombine
np.savetxt('noisecorr_{}_combine_M{}.txt'.format(args.thermostat,args.M), 
			np.column_stack((times, Kcombine*args.temperature, Kcombine/Kcombine[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time K K/K[0]')


# fsc=NoiseCorrSelfConsistent(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'], Kold=Kcombine, tstar=12)
