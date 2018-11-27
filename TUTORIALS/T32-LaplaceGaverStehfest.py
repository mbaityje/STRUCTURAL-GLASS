#!/bin/env python
################################################################
#
#
# DESCRIPTION
# This example shows how to Laplace transform and antitransform with the Gaver-Stehfest algorithm.
# The comments also describe some issues related to numerical precision: it is very easy to be in
# a parameter regime where nothing works because of numerical precision.
#
# To display help:
# python T32-LaplaceGaverStehfest.py -h
#
# To launch a simulation:
# python T32-LaplaceGaverStehfest.py -M5
#
# M is the Gaver-Stehfest parameter.
#
# Some references:
# http://www.unm.edu/~aierides/505/NotesOnNumericalLaplaceInversion.pdf
# https://www.kappaeng.com/PDF/Laplace_transform_numerical_inversion.pdf
# https://arxiv.org/pdf/1204.4754.pdf
# 
################################################################

#
# Propose several test functions, laplace transform them, and antitransform them back with the Gaver-Stehfest method
# 

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import scipy
from matplotlib import pyplot as plt
from lib.beylkin import Beylkin

import argparse
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-M','--M', type=int, required=False, default=5, help='number of fitting coefficients')
args = parser.parse_args()


def LaplaceSlow(f, p, tmax=np.inf):
	return quad(lambda t: f(t)*np.exp(-t*p), 0, tmax)[0]

def GaverStehfest(ftilde, t, M=args.M):
	if t==0:
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

#Test functions
def f1(t):
	return t*np.exp(-t)

def f1tilde(p):
	return 1./np.square(p+1)

left=0
right=10
n=100
delta=(right-left)/n
tarray=np.arange(left,right,delta)



# TEST LAPLACE TRANSFORM ON KNOWN RESULT
fig,ax=plt.subplots(1,1)
ax.set_xlim(left,right)
ax.plot(tarray, [LaplaceSlow(f1,xp) for xp in tarray], 'x', label='Laplace transform')
ax.plot(tarray, f1tilde(tarray), '+', label='Analytic solution')
ax.set_title('Laplace transform')
ax.set_xlabel('p')
ax.set_ylabel('$\\tilde{f}$ ($p$)')
ax.legend()
plt.show()

# TEST ANTI LAPLACE TRANSFORM ON KNOWN RESULT
fig,ax=plt.subplots(1,1)
ax.set_xlim(left,right)
ax.plot(tarray, [GaverStehfest(f1tilde,xp) for xp in tarray], '', label='Gaver-Stehfest')
ax.plot(tarray, f1(tarray), '+', label='Analytic solution')
ax.set_title('Laplace anti-transform')
ax.set_xlabel('t')
ax.set_ylabel('$f$ ($t$)')
ax.legend()
plt.show()



# NOW LET US WORK WITH ARRAYS INSTEAD OF FUNCTIONS
# Now tarray needs to be well-chosen so that the Gaver-Stehfest calculation
# does not ask for interpolations out of bounds.
# The maximum value of ftilde that needs to be calculated is
# t_max >= [k*log(2)/t]_max= 2*M*np.log(2)/tarray[1]
left=0
delta=tarray[1]
right=2*args.M*np.log(2)/delta
tarray=np.arange(left,right+delta,delta)

c1=f1(tarray)
c1tilde=f1tilde(tarray)


def chisquare(y,fx,delta):
	partials=np.square((y-fx)/delta)
	return partials.sum()

def FuncFromArray(x, y, kind='interp'):
	if kind == 'interp':
		return interp1d(x, y, kind='cubic')
	elif kind == 'fit':
		return FitFunction(x,y)

def FitFunction(x,y):
	ncoef=10
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function


# Now we create the interpolations. We have two ways to do so: through cubic spline interpolations, or through prony exponential fits.

interpc1=FuncFromArray(tarray, c1, kind='interp')
interpc1tilde=FuncFromArray(tarray, c1tilde, kind='interp')
fitc1=FuncFromArray(tarray, c1, kind='fit')
fitc1tilde=FuncFromArray(tarray, c1tilde, kind='fit')

# INTERPOLATION ERROR
# 
# As can be seen from the following figure, if the number of fitting coefficients is well-chosen 
# (ncoef not too large nor small -- in this example it is hardcoded),
# the approximation becomes much better than the interpolation (because we are interested in the maximum error, and because 
# everything below double numerical precision is the same for us).
# On the other side, the Gaver-Stehfest algorithm asks that the interpolation range be larger when M is larger. With exponentially
# decaying functions this can be a problem when we go fit them, since at some point the numerical noise takes over, and we are 
# fitting it (and it takes over the whole fitted function).
tarrayfine=np.arange(0,right,0.0025)
fig,ax=plt.subplots(1,1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(tarrayfine, np.abs(f1tilde(tarrayfine)-interpc1tilde(tarrayfine)), '', label='error on $\\tilde{f}_1(p)$ [interp]', color='darkblue')
ax.plot(tarrayfine, np.abs(f1tilde(tarrayfine)-fitc1tilde(tarrayfine)), '', label='error on $\\tilde{f}_1(p)$ [fit]', color='blue')
ax.plot(tarrayfine, np.abs(f1(tarrayfine)-interpc1(tarrayfine)), '', label='error on $f_1(t)$ [interp]', color='darkred')
ax.plot(tarrayfine, np.abs(f1(tarrayfine)-fitc1(tarrayfine)), '', label='error on $f_1(t)$ [fit]', color='red')
ax.set_title('Interpolation error')
ax.legend()
plt.show()



# TEST LAPLACE TRANSFORM ON KNOWN RESULT
# 
# Both fit and interpolation do well for the laplace transform
#
fig,ax=plt.subplots(1,1)
ax.set_xlim(left,right)
ax.plot(tarray, c1tilde, 'o', label='Analytic solution', color='grey')
ax.plot(tarray, [LaplaceSlow(interpc1, xp, tmax=tarray[-1]) for xp in tarray], 'x', label='Laplace transform [interp]', color='darkred')
ax.plot(tarray, [LaplaceSlow(fitc1, xp) for xp in tarray], '+', label='Laplace transform [fit]', color='red')
ax.set_title('Laplace transform')
ax.set_xlabel('p')
ax.set_ylabel('$\\tilde{f}$ ($p$)')
ax.legend()
plt.show()



# TEST GAVER-STEHFEST ANTI-TRANSFORM ON KNOWN RESULT
#
# The Gaver-Stefhest algorithm depends critically on the parameter M. M increases the precision, but only if the data is
# clean enough. As a rule of thumb, M can be at most the number of decimal digits of the error. Since the error on the interpolated
# function is around 1e-3, we can use maximum M=3, which is not very precise. With the fits we can get better precision, and go up to
# M=5,6,7, (best is 5) which is more or less the maximum M one can have with double precision (max is 8).
#
# If you replot the following figure in log-y scale you can see that also the fit does not do very well.
#
fig,ax=plt.subplots(1,1)
ax.set_xlim(left,right)
ax.plot(tarray, c1, 'o', label='Analytic solution', color='grey')
ax.plot(tarray, [GaverStehfest(interpc1tilde,xp) for xp in tarray], 'x', label='Gaver-Stehfest anti-transform [interp]', color='darkred')
ax.plot(tarray, [GaverStehfest(fitc1tilde,xp) for xp in tarray], '+', label='Gaver-Stehfest anti-transform [fit]', color='red')
ax.set_title('Laplace transform')
ax.set_xlabel('p')
ax.set_ylabel('$\\tilde{f}$ ($p$)')
ax.legend()
plt.show()


# I want to antitransform c1tilde and show that it is equal to c1
transform=np.array([LaplaceSlow(fitc1, xp) for xp in tarray])
fit_transform=FuncFromArray(tarray, transform, kind='fit')
antitransform=np.array([GaverStehfest(fit_transform, xp) for xp in tarray])

# TEST ANTI LAPLACE TRANSFORM ON KNOWN RESULT
fig,ax=plt.subplots(1,1)
ax.set_xlim(left,right)
ax.plot(tarray, antitransform, '', label='$L^{-1}[L[f(t)]]$')
ax.plot(tarray, c1, '+', label='$f(t)$')
ax.set_title('Laplace transform and anti-transform f1(x)')
ax.set_xlabel('t')
ax.set_ylabel('$f$ ($t$)')
ax.legend()
plt.show()


