#!/usr/bin/env python

import sys
import numpy as np
import argparse
import lib.module_measurements as med
from matplotlib import pyplot as plt
from lib.module_timelists import ListaLogaritmica
import lib.module_noisecorr as nc


# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('-M','--M', type=int, required=False, default=3, help='Number of particles')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--maxtime', type=float, required=False, default='-1.0', help='truncate the time axis at maxtime')
parser.add_argument('--tmin', type=float, required=False, default='0.05', help='time for the piecewise concatenation of Kvolterra with Klaplace')
parser.add_argument('--shiftCFP', action='store_true', help='if invoked, imposes CFP[0]=0')
parser.add_argument('--kind', required=False, choices=['interp','interp_lin','fit','combined'], default='combined', help='thermostat')
parser.add_argument('--ncoef', type=int, required=False, default=10, help='number of coefficients for fits')
parser.add_argument('--showplots', action='store_true', help='if invoked, shows plots during the calculation (this stops the calculation until the user closes the figure)')


args = parser.parse_args()

print(sys.argv)
print('T:',args.temperature)
print('thermostat:',args.thermostat)
print('kind:',args.kind)
print('M:',args.M)
print('ncoef:',args.ncoef)

if args.temperature<0:
	raise ValueError('Temperature='+str(args.temperature)+'. Must be positive. Aborting.')
else:
	invT=1./args.temperature



#READ KERNELS
CFF=np.load('CFF_'+str(args.thermostat)+'.npy') # access as CFP.item()['mean']
CFP=np.load('CFP_'+str(args.thermostat)+'.npy')
CPP=np.load('CPP_'+str(args.thermostat)+'.npy')
times=np.load('times_'+str(args.thermostat)+'.npy')

nt=len(times)
if not (len(CFP.item()['mean'])==nt and len(CFF.item()['mean'])==nt and len(CPP.item()['mean'])==nt):
	raise IndexError('The number of correlations ('+str(len(CFP.item()['mean']))+' or '+str(len(CFF.item()['mean']))+') is inconsistent with the number of times('+str(nt)+')')

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

#Put to zero CFP[0], to reduce a source of fluctuations
CFP.item()['mean'][0]-=CFP.item()['mean'][0] #Do not shift the whole curve, or there will be a bias in the integral!




#                         #
# HERE STARTS THE PROGRAM #
#                         #

cpp=np.copy(CPP.item()['mean'])
tlist=np.copy(times)

Klaplace=nc.NoiseCorrLaplace(times, cpp, args.M, args.temperature, ncoef=args.ncoef, kind=args.kind, showplots=args.showplots)
Kvolterra = nc.NoiseCorr(times=times, CFF=CFF.item()['mean'], CFP=CFP.item()['mean'], invT=invT)

if args.showplots:
        plt.title('Noise Correlation')
        plt.xlabel('$t$')
        plt.ylabel('$\mathcal{K}(t)$')
        plt.semilogx(times, Klaplace, label='$\\frac{kT-sC^P(s)}{C^P(s)}$')
        plt.semilogx(times, Kvolterra, label='K(t) [Volterra]')
        plt.legend()
        plt.show()

itmin=np.where(times>=args.tmin)[0][0]

Kcombine=np.copy(Kvolterra)
Kcombine[itmin:]=Klaplace[itmin:]


if args.showplots:
        plt.title('Noise Correlation')
        plt.semilogx(times, Kcombine, label='K(t) [Combine]')
        plt.legend()
        plt.show()

#Same Kcombine
np.savetxt('noisecorr_{}_combine_M{}.txt'.format(args.thermostat,args.M), 
			np.column_stack((times, Kcombine, Kcombine/Kcombine[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time K K/K[0]')
