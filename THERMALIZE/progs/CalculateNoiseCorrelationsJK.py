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
import lib.module_noisecorr as nc


# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
parser.add_argument('--thermostat', required=False, default='*', help='thermostat')
parser.add_argument('--maxtime', type=float, required=False, default='-1.0', help='truncate the time axis at maxtime')
parser.add_argument('-M','--M', type=int, required=False, default=3, help='2M = #of Gaver-Stehfest coefficients')
parser.add_argument('--showplots', action='store_true', help='If activated, shows a plot of computed correlations.')
parser.add_argument('--selfcon', action='store_true', help='If activated, calculates K also with the self-consistent method.')
parser.add_argument('--ncoef', type=int, required=False, default=10, help='number of coefficients for fits')
parser.add_argument('--tmin', type=float, required=False, default='0.05', help='time for the piecewise concatenation of Kvolterra with Klaplace')
parser.add_argument('--kind', required=False, choices=['interp','interp_lin','fit','combined'], default='combined', help='kind of interpolation')
parser.add_argument('--central_value', required=False, choices=['meanJK','medianJK','mean'], default='meanJK', help='what central value. meanJK (mean of JK blocks), medianJK (median of JK blocks), mean (K(t) calculated on the mean of CFF and CFP)')
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
	invT=np.float64(1./args.temperature)



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









#                         #
# HERE STARTS THE PROGRAM #
#                         #
K={'Volterra': {'blocksJK': np.ndarray((nblo,nt),dtype=np.float64)}, 'Laplace': {'blocksJK': np.ndarray((nblo,nt),dtype=np.float64)}, 'combine': {'blocksJK': np.ndarray((nblo,nt),dtype=np.float64)} }

for iblo in range(nblo):
	print('iblo = ',iblo)
	K['Laplace']['blocksJK'][iblo]=nc.NoiseCorrLaplace(times, CPP.item()['blocksJK'][iblo], args.M, args.temperature, ncoef=args.ncoef, kind=args.kind, showplots=args.showplots)
	K['Volterra']['blocksJK'][iblo]=nc.NoiseCorr(times, CFF.item()['blocksJK'][iblo], CFP.item()['blocksJK'][iblo], invT=invT, ncoef=args.ncoef)
	itmin, K['combine']['blocksJK'][iblo] = nc.NoiseCorrCombine(times, K['Volterra']['blocksJK'][iblo], K['Laplace']['blocksJK'][iblo], tmin=args.tmin)

if args.central_value=='medianJK':
	K['Laplace']['mean']=np.median(K['Laplace']['blocksJK'], axis=0)
	K['Volterra']['mean']=np.median(K['Volterra']['blocksJK'], axis=0)
elif args.central_value=='meanJK':
	K['Laplace']['mean']=np.mean(K['Laplace']['blocksJK'], axis=0)
	K['Volterra']['mean']=np.mean(K['Volterra']['blocksJK'], axis=0)
elif args.central_value=='mean':
	K['Laplace']['mean']=nc.NoiseCorrLaplace(times, CPP.item()['mean'], args.M, args.temperature, ncoef=args.ncoef, kind=args.kind, showplots=args.showplots)
	K['Volterra']['mean']=nc.NoiseCorr(times, CFF.item()['mean'], CFP.item()['mean'], invT=invT, ncoef=args.ncoef)


itmin, K['combine']['mean']=nc.NoiseCorrCombine(times, K['Volterra']['mean'], K['Laplace']['mean'])



K['Laplace' ]['errJK'] = np.sqrt(nblom1*(np.square(K['Laplace']['blocksJK']).mean(axis=0) - np.square(K['Laplace']['mean']) ) )
K['Volterra']['errJK'] = np.sqrt(nblom1*(np.square(K['Volterra']['blocksJK']).mean(axis=0) - np.square(K['Volterra']['mean'])))
K['combine']['errJK']=np.ndarray(nt,dtype=np.float64)
K['combine']['errJK'][:itmin]=K['Volterra']['errJK'][:itmin]
K['combine']['errJK'][itmin:]=K['Laplace' ]['errJK'][itmin:]




if args.showplots:
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







