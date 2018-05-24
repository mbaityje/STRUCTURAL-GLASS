#!/usr/bin/env python
################################################################
#
# Same as T19, but does not calculate CFF and CFP, in favor of reading them.
#
# DESCRIPTION
# This example:
# - reads CFF and CFP
# - calculates noise correlation function
# 
# Launch as:
# run T20-OrthogonalCorrelationsRead.py './test-output/correlations_dt0.0025_n10_T10_N1080_rho1.2.npz' --dt=0.0025 --temperature=10 --filter --truncate --scheme='rectangles' --label='prova'
# 
################################################################

from __future__ import print_function #for compatibility with python2
import sys #this has some system tools like sys.exit() that can be useful
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('corrname', help='file with CFF (.npz)')
parser.add_argument('-l','--label', required=False, default='', help='basename for the output files')
parser.add_argument('-T','--temperature', type=float, required=False, default=5, help='target Temperature')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='dt for MD integration')
parser.add_argument('--ncoef', type=int, required=False, default=0, help='number of fitting coefficients')
parser.add_argument('--filter', action='store_true', help='regularization that suppresses data with bad signal/noise ratio')
parser.add_argument('--truncate', action='store_true', help='set to zero the corr functions once the signal/noise ratio became too low')
parser.add_argument('--showfigs', action='store_true', help='whether to show figures (stops execution)')
parser.add_argument('--scheme', required=False, default='rectangles', help='integration scheme (retangles, trapeze, simpson)')
args = parser.parse_args()
print("corrname = ",args.corrname)
print("T = ",args.temperature)
print("dt = ",args.dt)
print("label = ",args.label)
print("ncoef = ",args.ncoef)
print("filter = ",args.filter)
print("truncate = ",args.truncate)
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
			break
	return myfit

def chisquare(y,fx,delta):
	partials=np.square((y-fx)/delta)
	return partials.sum()

from scipy.optimize import curve_fit
from lib.beylkinold import Beylkin
from lib.beylkin import Beylkin

corrfit={}
if args.ncoef==0:
	corrfit['CFF']=ManyFits(xdata,corr['CFF'],corr['errCFF'])
	corrfit['CFP']=ManyFits(xdata,corr['CFP'],corr['errCFP'])
else:
	bFF = Beylkin(decaying=True)
	bFF.driver_load(xdata, corr['CFF'], args.ncoef)
	corrfit['CFF']=np.append(bFF.correction(),bFF.correction()[-1])
	bFP = Beylkin(decaying=True)
	bFP.driver_load(xdata, corr['CFP'], args.ncoef)
	corrfit['CFP']=np.append(bFP.correction(),bFP.correction()[-1])


if args.showfigs:
	plt.xlim((0, 1))
	plt.errorbar(xdata[:len(corr   ['CFF'])], corr   ['CFF'], yerr=corr['errCFF'],errorevery=1, label='$\mathcal{C}^F_t/\mathcal{C}^F_0$', color='red')
	plt.errorbar(xdata[:len(corrfit['CFF'])], corrfit['CFF'],marker="o", label='$\mathcal{C}^F_t$ from fit points', color='blue')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t/\mathcal{C}_0$')
	plt.grid(True)
	plt.legend()
	plt.savefig('./test-output/corrFF'+args.label+".png")
	plt.show()

	plt.xlim((0, 1))
	plt.errorbar(xdata[:len(corr   ['CFP'])], corr['CFP'], yerr=corr['errCFP'],errorevery=1, label='$\mathcal{C}^F_t/\mathcal{C}^F_0$', color='red')
	plt.errorbar(xdata[:len(corrfit['CFP'])], corrfit['CFP'],marker="o", label='$\mathcal{C}^{FP}_t$ from fit', color='blue')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t/\mathcal{C}_0$')
	plt.grid(True)
	plt.legend()
	plt.savefig('./test-output/corrFP'+args.label+".png")
	plt.show()


################################################################
# FILTER THROUGH A REGULARIZATION FUNCTION
################################################################
def Truncate(func,err):
	'''
	Set to zero all the values of func after the signal became too small
	'''
	nsigma=1
	final=-1
	for i in range(len(func)):
		for j in range(i,len(func)):
			if func[j]>nsigma*err[j]:
				break
		if j==len(func)-1:
			final=i
			break
	if final>=0:
		print('Truncating function from the ',final,' element')
		func[final+1:]=0
	return


if args.showfigs:
	plt.errorbar(xdata[:len(corr   ['CFP'])], corr['CFP'], yerr=corr['errCFP'],errorevery=1, label='$\mathcal{C}^{FP}_t$ from data', color='green')
	plt.errorbar(xdata[:len(corrfit['CFP'])], corrfit['CFP'], label='$\mathcal{C}^{FP}_t$ from fit', color='lightgreen')

if args.filter:
	print("filter")
	corrfit['CFF']=np.tanh(np.abs(corrfit['CFF'])/corr['errCFF'])*corrfit['CFF']
	corrfit['CFP']=np.tanh(np.abs(corrfit['CFP'])/corr['errCFP'])*corrfit['CFP']
if args.truncate:
	print("Truncating")
	Truncate(corrfit['CFF'],corr['errCFF'])
	Truncate(corrfit['CFP'],corr['errCFP'])

if args.showfigs:
	plt.errorbar(xdata[:len(corrfit['CFP'])], corrfit['CFP'], label='$\mathcal{C}^{FP}_t$ regularized', color='red')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t$')
	plt.grid(True)
	plt.legend()
	plt.show()


################################################################
# WITH FITTED CORRELATION FUNCTIONS
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
maxiter=1000
Kold=np.copy(corr['CFF'])
Knew=np.zeros(Ncorr,dtype=np.float64)
xdata=np.arange(0,Ncorr)*args.dt

invT=np.float64(1./args.temperature)
for iter in range(maxiter):
	for n in range(Ncorr):
		Knew[n] = corrfit['CFF'][n] + invT * integral(corrfit['CFP'], Kold,n, args.dt,scheme=args.scheme)
	err=np.max(np.abs(Knew-Kold))
	print("iter:",iter," err = ",err)
	if err<1e-10: break
	Kold[:]=Knew[:] #Curiosamente, se metto Kold=Knew, converge subito al risultato giusto

plt.plot(xdata, Knew,label='$\mathcal{K}_t$')
plt.plot(xdata, corrfit['CFF'],label='$\mathcal{C}^{F}_t$')
plt.xlabel('$t$')
plt.ylabel('Correlation')
plt.grid(True)
plt.legend()
plt.title('Using fits - $T$ = %g, $dt$ = %g'%(args.temperature,args.dt))
plt.savefig('./test-output/corrNoise-selfconsistent'+args.label+"_KAfits.png")
plt.show()


