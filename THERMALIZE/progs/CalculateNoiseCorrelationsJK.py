#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
import lib.module_measurements as med
from matplotlib import pyplot as plt
from lib.beylkin import Beylkin
from scipy.integrate import quad
from lib.module_timelists import ListaLogaritmica



# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
# parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
# parser.add_argument('-L','--box_size', type=np.float64, required=True, help='Box Size (assumed cubic)')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--maxtime', type=float, required=False, default='-1.0', help='truncate the time axis at maxtime')
parser.add_argument('-M','--M', type=int, required=False, default=3, help='2M = #of Gaver-Stehfest coefficients')

args = parser.parse_args()

print(sys.argv)
print('M:',args.M)
# print('L:',args.box_size)
print('T:',args.temperature)
print('thermostat:',args.thermostat)
# L=args.box_size

# if not args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
# if not args.box_size>0: raise ValueError('L ({}) must be positive'.format(args.box_size))
if args.temperature<0:
	raise ValueError('Temperature='+str(args.temperature)+'. Must be positive. Aborting.')
else:
	invT=1./args.temperature


#READ KERNELS
CFF=np.load('CFFJK_'+str(args.thermostat)+'.npy') # access as CFP.item()['mean']
CFP=np.load('CFPJK_'+str(args.thermostat)+'.npy')
CPP=np.load('CPPJK_'+str(args.thermostat)+'.npy')
times=np.load('times_'+str(args.thermostat)+'.npy')
nt=len(times)
nblo=len(CFF.item()['blocks'])
nblom1=nblo-1
assert( nblo==len(CFF.item()['blocksJK']) and nblo==len(CFP.item()['blocks']) and nblo==len(CFP.item()['blocksJK']) )
it_one=np.where(times>1)[0][0]
if not (len(CFP.item()['mean'])==nt and len(CFF.item()['mean'])==nt):
	raise IndexError('The number of correlations ('+str(len(CFP.item()['mean']))+' or '+str(len(CFF.item()['mean']))+') is inconsistent with the number of times('+str(nt)+')')

invT=np.float64(1./args.temperature)

# PREPARE KERNELS
#Consider truncation of the time axis
if args.maxtime>0:
	times=np.array([t for t in times if t<args.maxtime])
	nt=len(times)
	CFF.item()['mean'    ]=CFF.item()['mean'    ][:nt]
	CFF.item()['errJK'   ]=CFF.item()['errJK'   ][:nt]
	CFF.item()['blocks'  ]=CFF.item()['blocks'  ][:nt]
	CFF.item()['blocksJK']=CFF.item()['blocksJK'][:nt]
	CFP.item()['mean'    ]=CFP.item()['mean'    ][:nt]
	CFP.item()['errJK'   ]=CFP.item()['errJK'   ][:nt]
	CFP.item()['blocks'  ]=CFP.item()['blocks'  ][:nt]
	CFP.item()['blocksJK']=CFP.item()['blocksJK'][:nt]
print('Integration interval: [{},{}]'.format(times[0],times[nt-1]))

#We know that CFP[0]=0
CFP.item()['mean'    ][0]    = 0
CFP.item()['errJK'   ][0]    = 0
for iblo in range(nblo):
	CFP.item()['blocks'  ][iblo][0] = 0
	CFP.item()['blocksJK'][iblo][0] = 0




def FitFunction(x, y, ncoef=10):
	'''Fit function as a sum of exponential functions'''
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function

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
		raise NotImplementedError("FuncFromArray: not implemented option {kind}")

def NoiseCorr(times, CFF, CFP, ncoef=12, kind='interp'):
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
	interpCFP= FuncFromArray(times, CFP, kind=kind, ncoef=ncoef)

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

def NoiseCorrBlocks(times, CFF, CFPblocks):
	'''
	Calculate noise correlation in the same way as NoiseCorr, but instead of interpolating on
	the average CFP, make the interpolation on every single block, and only later perform
	the average between the interpolations.
	The rationale is that perhaps the interpolation is introducing a bias that can be removed.
	In practice, the result doesn't change when this function is used, so there is no point in
	using this function.

	If you want to check this, just do

	KvolterraBlocks = NoiseCorrBlocks(times=times, CFF=CFF.item()['mean'], CFPblocks=CFP.item()['blocks'])

	'''
	nt=len(times)
	assert(nt==len(CFF) and nt==len(CFPblocks[0]))
	f=np.zeros(nt)
	g=CFF
	nblo=len(CFPblocks)

	interpCFP=[FuncFromArray(times, CFPblocks[iblo], kind='combined', ncoef=11) for iblo in range(nblo)]

	def kernel(i,j):
		if i==j: return 0
		dt=times[i]-times[j]
		av=np.array([interpCFP[iblo](dt) for iblo in range(nblo)]).mean()
		return av

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
	pmin=np.log(2)/x[-1]          -1e-8 #a small extra to make sure the range is enough in case of rounding errors 
	pmax=2*M*np.log(2)/(x[1]-x[0])+1e-8
	deltap=2*M*np.log(2)/(x[-1]-x[-2])
	print('pmin:',pmin,'pmax:',pmax, 'deltap:',deltap)
	p=np.arange(pmin, pmax+deltap, deltap)
	p=ListaLogaritmica(pmin, pmax+1e-3, 500)
	print('p[0 ] = ',p[0 ])
	print('p[-1] = ',p[-1])

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
		val=GaverStehfest(func_transform, xp, M)
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


def NoiseCorrLaplace(times, cpp, M, showplots=False, kind='combined'):
	pcpp,Lcpp,LLcpp=TransformAntitransform(times, cpp, M, ncoef=20, showplots=showplots, kind=kind)
	LK = np.array([( args.temperature - pcpp[ip]*Lcpp[ip] ) / Lcpp[ip] for ip in range(len(Lcpp))])

	LLK=np.ndarray(len(times))
	func_LK=FuncFromArray(pcpp, LK, kind='interp')

	for it in range(len(times)):
		tp=times[it]
		val=GaverStehfest(func_LK, tp, M)
		print('\rxp=%g,  LinvLf=%g'%(tp,val), end='')
		LLK[it]=val
	print('')
	return LLK

def NoiseCorrCombine(times, Kvolterra, Klaplace):
	'''
	Makes a step-wise K(t), equal to Kvolterra for short times and Klaplace for long times.
	'''
	tmin=0.02
	if args.temperature==5.0:
		tmin=0.027
	if args.temperature==2.0:
		tmin=0.018
	elif args.temperature==1.0:
		tmin=0.035
	elif args.temperature==0.8:
		tmin=0.04
	itmin=np.where(times>=tmin)[0][0]

	Kcombine=np.copy(Kvolterra)
	Kcombine[itmin:]=Klaplace[itmin:]
	return itmin, Kcombine

#                         #
# HERE STARTS THE PROGRAM #
#                         #
K={'Volterra': {'blocksJK': np.ndarray((nblo,nt),dtype=np.float64)}, 'Laplace': {'blocksJK': np.ndarray((nblo,nt),dtype=np.float64)}, 'combine': {} }

#
# Calculate K(t) by Laplace Transforming
# 
print('K(t) with the Laplace method')
for iblo in range(nblo):
	print('iblo = ',iblo)
	K['Laplace']['blocksJK'][iblo]=NoiseCorrLaplace(times, CPP.item()['blocksJK'][iblo], args.M, showplots=False, kind='combined')
K['Laplace']['mean']=K['Laplace']['blocksJK'].mean(axis=0)
K['Laplace']['errJK'] = np.sqrt(nblom1*(np.square(K['Laplace']['blocksJK']).mean(axis=0) - np.square(K['Laplace']['mean']) ) )

#
# Calculate K(t) by integrating the Volterra Equation
# 
print('K(t) with the Volterra method')
ncoef=12
for iblo in range(nblo):
	print('iblo = ',iblo)
	K['Volterra']['blocksJK'][iblo]=NoiseCorr(times, CFF.item()['blocksJK'][iblo], CFP.item()['blocksJK'][iblo], ncoef=ncoef)/args.temperature
K['Volterra']['mean']=K['Volterra']['blocksJK'].mean(axis=0)
K['Volterra']['errJK'] = np.sqrt(nblom1*(np.square(K['Volterra']['blocksJK']).mean(axis=0) - np.square(K['Volterra']['mean']))) #Togliere np.abs


#
# Concatenate Kvolterra and Klaplace
# 
print('Concatenating K(t)')
itmin, K['combine']['mean']=NoiseCorrCombine(times, K['Volterra']['mean'], K['Laplace']['mean'])
K['combine']['blocksJK']=np.ndarray((nblo,nt), dtype=np.float64)
for iblo in range(nblo):
	itmin, K['combine']['blocksJK'][iblo] = NoiseCorrCombine(times, K['Volterra']['blocksJK'][iblo], K['Laplace']['blocksJK'][iblo])

K['combine']['errJK']=np.ndarray(nt,dtype=np.float64)
K['combine']['errJK'][:itmin]=K['Volterra']['errJK'][:itmin]
K['combine']['errJK'][itmin:]=K['Laplace' ]['errJK'][itmin:]


fig,ax=plt.subplots(1,1)
ax.set_xscale('log')
ax.errorbar(times, K['Volterra']['mean'], yerr=K['Volterra']['errJK'], label='Volterra')
ax.errorbar(times, K['Laplace' ]['mean'], yerr=K['Laplace' ]['errJK'], label='Laplace')
ax.grid(axis='y', color='black')
ax.legend()
plt.show()

np.savetxt('noisecorrJK_{}_M{}.txt'.format(args.thermostat,args.M), 
			np.column_stack((
				times, 
				K['combine']['mean'], K['combine']['errJK'], 
				K['combine']['mean']/K['combine']['mean'][0], 
				np.abs(K['combine']['errJK']/K['combine']['mean'][0])+ K['combine']['errJK'][0]*K['combine']['mean']/np.square(K['combine']['mean'][0]),
				K['Volterra']['mean'], K['Volterra']['errJK'], 
				K['Volterra']['mean']/K['Volterra']['mean'][0], 
				np.abs(K['Volterra']['errJK']/K['Volterra']['mean'][0])+ K['Volterra']['errJK'][0]*K['Volterra']['mean']/np.square(K['Volterra']['mean'][0]),
				K['Laplace']['mean'], K['Laplace']['errJK'], 
				K['Laplace']['mean']/K['Laplace']['mean'][0], 
				np.abs(K['Laplace']['errJK']/K['Volterra']['mean'][0])+ K['Volterra']['errJK'][0]*K['Laplace']['mean']/np.square(K['Volterra']['mean'][0])
				)), 
			fmt=['%g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g','%.14g'], 
			header='1)time 2-3)Kcomb 4-5)Kcomb/Kcomb[0] 6-7)Kvolt 8-9)Kvolt/Kvolt[0] 10-11)Klaplace 12-13)Klaplace/Kvolt[0]'
			)

np.save('noisecorrJK_{}_M{}.npy'.format(args.thermostat,args.M), K)



