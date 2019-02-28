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
parser.add_argument('corrname', help='correlation function we want to read')
parser.add_argument('--temperature', type=float, required=True, help='temperature is needed to transform the correlator into friction')
parser.add_argument('--thermostat', required=True, help='thermostat is needed for output files')
parser.add_argument('--density', type=float, required=False, default=1.2, help='density is optional because I always work at rho=1.2')
parser.add_argument('--tstar', type=float, required=False, default=None, help='no integration for t>tstar')
parser.add_argument('--showplots', action='store_true', help='If activated, plots are shown')
args = parser.parse_args()

if (not args.tstar==None) and args.tstar<=0: raise ValueError('tstar should be positive ({})'.format(args.tstar))

data=np.loadtxt(args.corrname)
times=data[:,0]
C=data[:,1]
# errC=data[:,2]
del data
istar= -1 if args.tstar==None else np.where(times>args.tstar)[0][0]




#SHORT-TIME FUNCTION
ishort=np.where(times>0.01)[0][0]
def fshort(x, a1, a2):
	return a1/np.cosh(a2*x)
params=[C[0],1]
params,dummy=curve_fit(fshort, times[:ishort], C[:ishort], p0=params)
if args.showplots:
	plt.subplot(111, xscale="log", yscale="linear", ylim=(0,1+C[0]))
	plt.plot(times, fshort(times,params[0],params[1]), times, C)
	plt.show()


#LONG-TIME FUNCTION
Clong=C-fshort(times,params[0], params[1])
#A little more cleaning due to the fact that we attached two functions
Clong[np.where(times<=0.05)]=0
Clong[np.where(Clong<0)]=0


if args.showplots:
	plt.subplot(111, xscale="log", yscale="linear")
	plt.plot(times, Clong, times, C)
	plt.show()


#FRICTION
total_friction = np.trapz(C[:istar], x=times[:istar]) * args.density/args.temperature
friction_short=quad(fshort, 0, 1, args=(params[0],params[1]))[0] #At t=1 the short time regime is over
friction_long=np.trapz(Clong[:istar], x=times[:istar]) * args.density/args.temperature


##############################################
# 
# OUTPUT
#
outdir=path.dirname(args.corrname)+'/'
print('outdir = ',outdir)

#
# Friction
#
out_friction = open(outdir+'frictionNoise_'+args.thermostat+'.txt',"w") 
out_friction.write('#total_friction != friction_long+friction_short since because of the piecewise attachment I set friction_long to zero for short times.\n')
out_friction.write('#total_friction friction_long friction_short\n')
out_friction.write('{} {} {} {} {}\n'.format(args.temperature,args.thermostat,total_friction,friction_long,friction_short))
out_friction.close()

#
# Fit parameters
#
out_short = open(outdir+'shortNoise_'+args.thermostat+'.txt',"w")
out_short.write("# Fit parameters for short-time behavior. fshort(x)=a1/cosh(a2)\n")
out_short.write("#T thermostat a1 a2\n")
out_short.write("{} {} {} {}\n".format(args.temperature,args.thermostat,params[0],params[1]))
out_short.close()

#
# Long-time correlation function
#
base,extension=path.splitext(path.basename(args.corrname))
np.savetxt(base+'_long'+extension, 
			np.column_stack((times, Clong, Clong/C[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time Klong Klong/K[0]')


