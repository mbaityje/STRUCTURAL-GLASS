#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import curve_fit
import lib.module_measurements as med
from matplotlib import pyplot as plt
from os import path

# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('nametxt', help='correlation function we want to read')
parser.add_argument('nameblocks', help='JK blocks of the correlation function we want to read')
parser.add_argument('--temperature', type=float, required=True, help='temperature is needed to transform the correlator into friction')
parser.add_argument('--thermostat', required=True, help='thermostat is needed for output files')
parser.add_argument('--density', type=float, required=False, default=1.2, help='density is optional because I always work at rho=1.2')
parser.add_argument('--tstar', type=float, required=False, default=None, help='no integration for t>tstar')
parser.add_argument('--showplots', action='store_true', help='If activated, plots are shown')
args = parser.parse_args()

if (not args.tstar==None) and args.tstar<=0: raise ValueError('tstar should be positive ({})'.format(args.tstar))

K=np.load(args.nameblocks).item()
data=np.loadtxt(args.nametxt)
times=data[:,0]
C=data[:,1]
errC=data[:,2]
istar= -1 if args.tstar==None else np.where(times>args.tstar)[0][0]


nt=len(times)
nblo=len(K['combine']['blocksJK'])
nblom1=nblo-1

obs={
	'total_friction': {'blocksJK': np.ndarray(nblo,dtype=np.float64)}, 
	'friction_short': {'blocksJK': np.ndarray(nblo,dtype=np.float64)}, 
	'friction_long' : {'blocksJK': np.ndarray(nblo,dtype=np.float64)},
	'shortParams'  : {'blocksJK': np.ndarray((nblo,2),dtype=np.float64)},
	}


#SHORT-TIME FUNCTION PREP
ishort=np.where(times>0.01)[0][0]
def fshort(x, a1, a2):
	return a1/np.cosh(a2*x)


#
for iblo in range(nblo):

	#SHORT-TIME FUNCTION
	C=K['combine']['blocksJK'][iblo]
	params=[C[0],1]
	params,dummy=curve_fit(fshort, times[:ishort], C[:ishort], p0=params)
	obs['shortParams']['blocksJK'][iblo]=params


	#LONG-TIME FUNCTION
	Clong=C-fshort(times,params[0], params[1])
	#A little more cleaning due to the fact that we attached two functions
	Clong[np.where(times<=0.05)]=0
	Clong[np.where(Clong<0)]=0

	#FRICTION
	obs['total_friction']['blocksJK'][iblo] = np.trapz(C[:istar], x=times[:istar]) * args.density/args.temperature
	obs['friction_short']['blocksJK'][iblo] = quad(fshort, 0, 1, args=(params[0],params[1]))[0] #At t=1 the short time regime is over
	obs['friction_long' ]['blocksJK'][iblo] = np.trapz(Clong[:istar], x=times[:istar]) * args.density/args.temperature

	if args.showplots:
		plt.subplot(111, xscale="log", yscale="linear")
		plt.plot(times, C)

if args.showplots: plt.show()


#Mean
C=K['combine']['mean']
params=[C[0],1]
params,dummy=curve_fit(fshort, times[:ishort], C[:ishort], p0=params)
if args.showplots:
	plt.subplot(111, xscale="log", yscale="linear", ylim=(0,1+C[0]))
	plt.plot(times, fshort(times,params[0],params[1]), times, C)
	plt.show()
obs['shortParams']['mean']=params

Clong=C-fshort(times,params[0], params[1])
#A little more cleaning due to the fact that we attached two functions
Clong[np.where(times<=0.05)]=0
Clong[np.where(Clong<0)]=0
if args.showplots:
	plt.subplot(111, xscale="log", yscale="linear")
	plt.plot(times, Clong, times, C)
	plt.show()

obs['total_friction']['mean'] = np.trapz(C[:istar], x=times[:istar]) * args.density/args.temperature
obs['friction_short']['mean'] = quad(fshort, 0, 1, args=(params[0],params[1]))[0] #At t=1 the short time regime is over
obs['friction_long' ]['mean'] = np.trapz(Clong[:istar], x=times[:istar]) * args.density/args.temperature

#Error
obs['shortParams']['errJK'] = np.sqrt(nblom1*(np.square(obs['shortParams']['blocksJK']).mean(axis=0) - np.square(obs['shortParams']['mean']) ) )

obs['total_friction']['errJK'] = np.sqrt(nblom1*(np.square(obs['total_friction']['blocksJK']).mean(axis=0) - np.square(obs['total_friction']['mean']) ) )
obs['friction_short']['errJK'] = np.sqrt(nblom1*(np.square(obs['friction_short']['blocksJK']).mean(axis=0) - np.square(obs['friction_short']['mean']) ) )
obs['friction_long' ]['errJK'] = np.sqrt(nblom1*(np.square(obs['friction_long' ]['blocksJK']).mean(axis=0) - np.square(obs['friction_long' ]['mean']) ) )


##############################################
# 
# OUTPUT
#
outdir=path.dirname(args.nametxt)+'/'
print('outdir = ',outdir)

#
# Friction
#
out_friction = open(outdir+'frictionNoiseJK_'+args.thermostat+'.txt',"w") 
out_friction.write('#T thermostat total_friction err friction_long err friction_short err\n')
out_friction.write('{} {} {} {} {} {} {} {}\n'.format(args.temperature,args.thermostat,
	obs['total_friction']['mean'],obs['total_friction']['errJK'],
	obs['friction_long' ]['mean'],obs['friction_long' ]['errJK'],
	obs['friction_short']['mean'],obs['friction_short']['errJK']))
out_friction.close()


#
# Fit parameters
#
out_short = open(outdir+'shortNoiseJK_'+args.thermostat+'.txt',"w")
out_short.write("# Fit parameters for short-time behavior. fshort(x)=a1/cosh(a2)\n")
out_short.write("#T thermostat a1 err a2 err\n")
out_short.write("{} {} {} {} {} {}\n".format(args.temperature,args.thermostat,
	obs['shortParams']['mean'][0],obs['shortParams']['errJK'][0],
	obs['shortParams']['mean'][1],obs['shortParams']['errJK'][1]))
out_short.close()

