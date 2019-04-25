#!/usr/bin/env python
#
# Module for time lists
# 
from __future__ import print_function
import numpy as np
from lib.module_timelists import ListaLogaritmica
from scipy.interpolate import interp1d
from lib.beylkin import Beylkin
from matplotlib import pyplot as plt
from scipy.integrate import quad

# FITTING AND INTERPOLATION

def FitFunction(x, y, ncoef=10):
	'''Fit function as a sum of exponential functions'''
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function

def FuncFromArray(x, y, kind='interp', ncoef=10):
	'''
	Several ways of obtaining a smooth function (needed e.g. to calculate the laplace transforms
	and antitransforms) from an array (i.e. the functions that we have as input).

	interp: cubic spline interpolation.
	
	interp_lin: linear spline interpolation.
	
	fit: a fit with the Prony scheme. Since Prony is implemented only on a linear grid, be 
	careful of applying it only when the grid of times is linear.

	combined: this is my dummy way to deal with a non-linear grid (my data is on 
	non-linear grids). I first do a cubic spline on the non-linear grid, use the spline to
	make it linear, and then use Prony on the resulting linear grid. It is not too
	slow and the fits are really good.
	'''
	if kind == 'interp':
		return interp1d(x, y, kind='cubic', assume_sorted=True)
	elif kind == 'interp_lin':
		return interp1d(x, y, kind='slinear', assume_sorted=True)
	elif kind == 'fit':
		return FitFunction(x,y, ncoef=ncoef)
	elif kind == 'combined':
		tempfunc=interp1d(x, y, kind='cubic', assume_sorted=True)
		it_one = np.where(x>=1)[0][0]
		mydelta=(x[it_one]-x[0])/1000
		xlinear=np.arange(x[0], x[it_one]+mydelta, mydelta)
		return FitFunction(xlinear,tempfunc(xlinear), ncoef=ncoef)
	else:
		raise NotImplementedError('FuncFromArray: not implemented option')

###########
# LAPLACE #
###########


def LaplaceSlow(f, p, tmax=np.inf):
	return quad(lambda t: f(t)*np.exp(-t*p), 0, tmax)[0]


def GaverStehfest(ftilde, t, M):
	''' Inverse Laplace transform with the Gaver-Stehfest method	'''
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



######################
# NOISE CORRELATIONS #
######################
#
# K = <F(0)F(t*)/T>, where t* is evolved through the orthogonal dynamics.
#


def NoiseCorr(times, CFF, CFP, invT, ncoef=12, kind='interp'):
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
			print('NoiseCorr is giving unphysical values (it is growing)')
	print('')
	return f*invT


def NoiseCorrBlocks(times, CFF, CFPblocks, invT):
	'''
	Calculate noise correlation in the same way as NoiseCorr, but instead of interpolating on
	the average CFP, make the interpolation on every single block, and only later perform
	the average between the interpolations.
	The rationale is that perhaps the interpolation is introducing a bias that can be removed.
	In practice, the result doesn't change when this function is used, so there is no point in
	using this function.

	If you want to check this, just do

	KvolterraBlocks = NoiseCorrBlocks(times=times, CFF=CFF.item()['mean'], CFPblocks=CFP.item()['blocks'], invT=1./args.temperature)

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
	return f*invT


def NoiseCorrLaplace(times, cpp, M, T, ncoef=10, kind='combined', showplots=False):
	pcpp,Lcpp,LLcpp=TransformAntitransform(times, cpp, M, ncoef=ncoef, showplots=showplots, kind=kind)
	LK = np.array([( T - pcpp[ip]*Lcpp[ip] ) / Lcpp[ip] for ip in range(len(Lcpp))])

	LLK=np.ndarray(len(times))
	func_LK=FuncFromArray(pcpp, LK, kind='interp')

	for it in range(len(times)):
		tp=times[it]
		val=GaverStehfest(func_LK, tp, M)
		print('\rxp=%g,  LinvLf=%g'%(tp,val), end='')
		LLK[it]=val
	print('')
	return LLK


def NoiseCorrCombine(times, Kvolterra, Klaplace, tmin=0.05):
	'''
	Makes a step-wise K(t), equal to Kvolterra for short times and Klaplace for long times.
	'''
	itmin=np.where(times>=tmin)[0][0]
	Kcombine=np.copy(Kvolterra)
	Kcombine[itmin:]=Klaplace[itmin:]
	return itmin, Kcombine

def NoiseCorrSelfConsistent(times, CFF, CFP, invT, ncoef=10, Kold=None, maxiter=1000, tstar=None, fits=False):
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

	interpCFP  = FitFunction(times, CFP , ncoef=ncoef) if fits else interp1d(times, CFP , kind='cubic', assume_sorted=True)
	interpKold = FitFunction(times, Kold, ncoef=ncoef) if fits else interp1d(times, Kold, kind='cubic', assume_sorted=True)

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
		interpKold=FitFunction(times, Knew, ncoef=ncoef) if fits else interp1d(times, Knew, kind='cubic', assume_sorted=True)
	return Knew*invT


if __name__=='__main__':
	pass