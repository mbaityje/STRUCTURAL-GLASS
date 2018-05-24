#!/usr/bin/env python
################################################################
#
# Same as T21, instead of fits I use interpolations.
# This allows strightforwardly for extrapolations and interpolations.
# Also, the integration grid can be made more or less dense throught --fact argument.
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
parser.add_argument('--dt'  , type=float, required=False, default=0.0025, help='dt for MD integration')
parser.add_argument('--fact', type=float, required=False, default=1.0   , help='if >1, the integration grid is more dense, if <1 the grid is less dense. [default:1]')
parser.add_argument('--showfigs', action='store_true', help='whether to show figures (stops execution)')
parser.add_argument('--scheme', required=False, default='rectangles', help='integration scheme (retangles, trapeze, simpson)')
args = parser.parse_args()
print("corrname = ",args.corrname)
print("T = ",args.temperature)
print("dt = ",args.dt)
print("label = ",args.label)
print("scheme = ",args.scheme)
print("fact = ",args.fact)
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
dtfine=args.dt/args.fact
invT=np.float64(1./args.temperature)
xmax=xdata[-1]
xdatafine=np.arange(0, xmax, args.dt/args.fact)
NcorrFine=len(xdatafine)
du=xdatafine[1]-xdatafine[0]


################################################################
# FIT CORRELATION FUNCTIONS 
################################################################
from scipy.interpolate import interp1d
from scipy.integrate import quad

interpCFF=interp1d(xdata, corr['CFF'], kind='cubic')
interpCFP=interp1d(xdata, corr['CFP'], kind='cubic')


if args.showfigs:
	# plt.xlim((0, 1))
	plt.errorbar(xdata, corr   ['CFF'], yerr=corr['errCFF'],errorevery=1, label='$\mathcal{C}^F_t$', color='red')
	plt.errorbar(xdatafine, interpCFF(xdatafine), label='$\mathcal{C}^F_t$ interpolated', color='green')
	plt.xlabel('$t$')
	plt.ylabel('$\mathcal{C}_t$')
	plt.grid(True)
	plt.legend()
	plt.savefig('./test-output/corrFF'+args.label+".png")
	plt.show()

	# plt.xlim((0, 1))
	plt.errorbar(xdata, corr['CFP'], yerr=corr['errCFP'],errorevery=1, label='$\mathcal{C}^{FP}_t$', color='red')
	plt.errorbar(xdatafine, interpCFP(xdatafine), label='$\mathcal{C}^{FP}_t$ interpolated', color='green')
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
print('Use a different dt than the one of the measurements, du=',du)
CFP=np.copy(interpCFP(xdatafine))-corr['CFP'][0]
CFF=np.copy(interpCFF(xdatafine))#-corr['CFF'][200:].mean()
Kold=np.copy(interpCFF(xdatafine))#-corr['CFF'][200:].mean()
Knew=np.zeros(NcorrFine,dtype=np.float64)

maxiter=1000
for iter in range(maxiter):
	for n in range(NcorrFine):
		Knew[n] = CFF[n] + invT * integral(CFP, Kold,n, dtfine,scheme=args.scheme)
	err=np.max(np.abs(Knew-Kold))
	print("iter:",iter," err = ",err)
	if err<1e-10: break
	Kold[:]=Knew[:]

plt.plot(xdatafine, Knew,label='$\mathcal{K}_t$')
plt.plot(xdatafine, interpCFF(xdatafine),label='$\mathcal{C}^{F}_t$')
plt.xlabel('$t$')
plt.ylabel('Correlation')
plt.grid(True)
plt.legend()
plt.title('Using fits - $T$ = %g, $dt$ = %g, $du$ = %g'%(args.temperature,args.dt,du))
plt.savefig('./test-output/corrNoise-selfconsistent'+args.label+"_KAfits.png")
plt.show()

