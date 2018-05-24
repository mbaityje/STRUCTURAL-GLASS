#!/usr/bin/env python
################################################################
#
# Same as T20, but the fits are represented through a function instead of an array.
# This allows strightforwardly for extrapolations and interpolations.
#
# DESCRIPTION
# This example:
# - reads CFF and CFP
# - performs fits
# - calculates noise correlation function
# 
# Launch as:
# run T21-OrthogonalCorrelationsFancyFits.py './test-output/correlations_dt0.0025_n10_T10_N1080_rho1.2.npz' --dt=0.0025 --temperature=10 --filter --truncate --scheme='rectangles' --label='prova'
# 
################################################################

from __future__ import print_function #for compatibility with python2
import sys #this has some system tools like sys.exit() that can be useful
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
import matplotlib.pyplot as plt
from lib.beylkin import Beylkin
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('corrname', help='file with CFF (.npz)')
parser.add_argument('-l','--label', required=False, default='', help='basename for the output files')
parser.add_argument('-T','--temperature', type=float, required=False, default=5, help='target Temperature')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='dt for MD integration')
parser.add_argument('--ncoef', type=int, required=False, default=0, help='number of fitting coefficients')
parser.add_argument('--crop', type=int, required=False, default=0, help='how many points to use for the fits (0: Ncorr [default])')
parser.add_argument('--showfigs', action='store_true', help='whether to show figures (stops execution)')
parser.add_argument('--scheme', required=False, default='rectangles', help='integration scheme (retangles, trapeze, simpson)')
args = parser.parse_args()
print("corrname = ",args.corrname)
print("T = ",args.temperature)
print("dt = ",args.dt)
print("label = ",args.label)
print("ncoef = ",args.ncoef)
print("crop = ",args.crop)
print("scheme = ",args.scheme)
print("showfigs = ",args.showfigs)
assert(args.temperature>0)

################################################################
# READ STANDARD CORRELATORS
################################################################

#Save correlation functions
corr=np.load(args.corrname)
Ncorr=len(corr['CFF'])
if len(corr['errCFF'])!= Ncorr: raise ValueError('Inconsistent data')
if len(corr[   'CFP'])!= Ncorr: raise ValueError('Inconsistent data')
if len(corr['errCFP'])!= Ncorr: raise ValueError('Inconsistent data')
xdata=np.arange(0,Ncorr)*args.dt

#The finer grid
factor=10
dtfine=args.dt/factor
NcorrFine=Ncorr*factor
invT=np.float64(1./args.temperature)
xdatafine=np.arange(0, Ncorr*args.dt, args.dt/factor)
crop=Ncorr if args.crop==0 else args.crop #Make fit over the first part of the curve


################################################################
# FIT CORRELATION FUNCTIONS 
################################################################
def ManyFits(xdata,function,err, maxncoef=50, crop=0):
	'''Fit with increasing number of parameters, until chi^2/ndof<1'''
	max_x=Ncorr if crop==0 else crop
	for ncoef in range(1,maxncoef):
		ndof=max_x-ncoef-1
		b = Beylkin(decaying=True)
		b.driver_load(xdata[0:max_x], function[0:max_x], ncoef)
		myfit=np.append(b.correction(),b.correction()[-1])
		chisq=chisquare(function[0:max_x], myfit, err[0:max_x])
		print('n: %d,  chi^2 = %g'%(ncoef,chisq))
		if chisq<ndof:
			if np.max(b.prony_function(xdatafine))<2*np.max(function): #A loose check to prevent overfitting
				break
	return myfit,b.prony_function

def chisquare(y,fx,delta):
	partials=np.square((y-fx)/delta)
	return partials.sum()

corrfit={}
if args.ncoef==0:
	corrfit['CFF'],funcCFF=ManyFits(xdata,corr['CFF'],corr['errCFF'],crop=crop)
	corrfit['CFP'],funcCFP=ManyFits(xdata,corr['CFP'],corr['errCFP'],crop=crop)
else:
	bFF = Beylkin(decaying=True)
	bFF.driver_load(xdata, corr['CFF'], args.ncoef)
	corrfit['CFF']=np.append(bFF.correction(),bFF.correction()[-1])
	funcCFF=bFF.prony_function
	bFP = Beylkin(decaying=True)
	bFP.driver_load(xdata, corr['CFP'], args.ncoef)
	corrfit['CFP']=np.append(bFP.correction(),bFP.correction()[-1])
	funcCFP=bFP.prony_function


if args.showfigs:
	plt.xlim((0, 1))
	plt.errorbar(xdata, corr   ['CFF'], yerr=corr['errCFF'],errorevery=1, label='$\mathcal{C}^F_t$', color='red')
	plt.errorbar(xdata[:len(corrfit['CFF'])], corrfit['CFF'],marker="o", label='$\mathcal{C}^F_t$ from fit points', color='blue')
	plt.errorbar(xdatafine, funcCFF(xdatafine), label='$\mathcal{C}^F_t$ from fit function', color='green')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t$')
	plt.grid(True)
	plt.legend()
	plt.savefig('./test-output/corrFF'+args.label+".png")
	plt.show()

	plt.xlim((0, 1))
	plt.errorbar(xdata, corr['CFP'], yerr=corr['errCFP'],errorevery=1, label='$\mathcal{C}^{FP}_t$', color='red')
	plt.errorbar(xdata[:len(corrfit['CFP'])], corrfit['CFP'],marker="o", label='$\mathcal{C}^{FP}_t$ from fit', color='blue')
	plt.errorbar(xdatafine, funcCFP(xdatafine), label='$\mathcal{C}^{FP}_t$ from fit function', color='green')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t$')
	plt.grid(True)
	plt.legend()
	plt.savefig('./test-output/corrFP'+args.label+".png")
	plt.show()



################################################################
# 
# MEASURE K WITH FITTED CORRELATION FUNCTIONS
# 
# If I integrate things don't go well because the dt of the measurements is large.
# To make the dt thinner, we 
# - define a new thinner grid
# 
# 
################################################################
from scipy.integrate import simps,romb
def integral(myCFP, myKold, myn, dt, scheme='rectangles'):
	'''Integral for the self-consistent calculation'''
	if n==0: return 0
	integrand=np.array([myCFP[myn-i]*myKold[i] for i in range(myn)])

	if scheme   == 'rectangles':
		output=integrand.sum()*dt
	elif scheme == 'trapeze':
		output=np.trapz(integrand,dx=dt)
	elif scheme == 'simpson':
		output=simps(integrand,dx=dt, even='avg')
	else:
		raise NotImplementedError('Integration scheme must be rectangles, trapeze or simpson. '+scheme+' is not a valid option.')
	return output



print('Measure Noise correlation functions')
print('integrate with the '+args.scheme+' method')
print('Use a finer dt than the one of the measurements')
Kold=np.copy(funcCFF(xdatafine))
Knew=np.zeros(NcorrFine,dtype=np.float64)
CFP=np.copy(funcCFP(xdatafine))
CFF=np.copy(funcCFF(xdatafine))

maxiter=1000
for iter in range(maxiter):
	for n in range(NcorrFine):
		Knew[n] = CFF[n] + invT * integral(CFP, Kold,n, dtfine,scheme=args.scheme)
	err=np.max(np.abs(Knew-Kold))
	print("iter:",iter," err = ",err)
	if err<1e-10: break
	Kold[:]=Knew[:]

plt.plot(xdatafine, Knew,label='$\mathcal{K}_t$')
plt.plot(xdatafine, funcCFF(xdatafine),label='$\mathcal{C}^{F}_t$')
plt.xlabel('$t$')
plt.ylabel('Correlation')
plt.grid(True)
plt.legend()
plt.title('Using fits - $T$ = %g, $dt$ = %g'%(args.temperature,args.dt))
plt.savefig('./test-output/corrNoise-selfconsistent'+args.label+"_KAfits.png")
plt.show()


