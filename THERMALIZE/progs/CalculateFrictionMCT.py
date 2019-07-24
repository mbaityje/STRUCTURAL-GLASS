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
parser.add_argument('Kname', help='correlation function we want to read')
parser.add_argument('KMCTname', help='correlation function we want to read')
parser.add_argument('--temperature', type=float, required=True, help='temperature is needed to transform the correlator into friction')
parser.add_argument('--thermostat', required=True, help='thermostat is needed for output files')
parser.add_argument('--density', type=float, required=False, default=1.2, help='density is optional because I always work at rho=1.2')
parser.add_argument('--showplots', action='store_true', help='If activated, plots are shown')
args = parser.parse_args()

#Read K(t)
data=np.loadtxt(args.Kname)
times=data[:,0]
K=data[:,1]
Kerr=data[:,2]

#Read Kmct(t)
data=np.loadtxt(args.KMCTname)
timesMCT=data[:,0]
kMCT=data[:,1]
kMCTerr=data[:,2]

del data
istar= -1


def fshort(x, a1, a2):
	return a1/np.cosh(a2*x)

def FitShortTime(times, C, endShort=0.1):
	'''Fits K(t) or Kmct(t) at short times'''

	ishort=np.where(times>endShort)[0][0]

	params=[C[0],1]
	params,cov=curve_fit(fshort, times[:ishort], C[:ishort], p0=params)
	if args.showplots:
		plt.subplot(111, xscale="log", yscale="linear", ylim=(0,1+C[0]))
		plt.plot(times, fshort(times,params[0],params[1]), label='fit')
		plt.plot(times, C, label='data')
		plt.legend()
		plt.show()

	#Long-time function by subtracting the short time
	Clong=C-fshort(times,params[0], params[1])
	#A little more cleaning due to the fact that we attached two functions
	ini_long=np.where(Clong[:int(len(Clong)/3)]<0)[0][-1]
	Clong[0:ini_long+1]=0

	print('times:',times)

	if args.showplots:
		plt.subplot(111, xscale="log", yscale="linear")
		plt.plot(times, Clong, times, C)
		plt.show()
	return params, cov, times, Clong

tempK   =FitShortTime(times, K)
tempKmct=FitShortTime(timesMCT, kMCT)
fit={'K'   :{'params': tempK   [0], 'cov': tempK   [1], 'times': tempK   [2], 'Clong': tempK   [3]} ,
	 'Kmct':{'params': tempKmct[0], 'cov': tempKmct[1], 'times': tempKmct[2], 'Clong': tempKmct[3]}
	 }
del tempK, tempKmct


# Plot long time of Kmct and short time of the exact
plt.subplot(111, xscale="log", yscale="linear")
plt.plot(fit['K']['times'], fshort(fit['K']['times'], fit['K']['params'][0], fit['K']['params'][1]), label='Kshort')
plt.plot(fit['Kmct']['times'], fit['Kmct']['Clong'], label='KmctLong')
plt.show()


sys.exit()




#LONG-TIME FUNCTION




#FRICTION
total_friction = np.trapz(C[:istar], x=times[:istar]) * args.density/args.temperature
friction_short=quad(fshort, 0, 1, args=(params[0],params[1]))[0] #At t=1 the short time regime is over
friction_long=np.trapz(KmctLong[:istar], x=times[:istar]) * args.density/args.temperature


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
out_short.write("#T thermostat a1 err1 a2 err2\n")
out_short.write("{} {} {} {} {} {}\n".format(args.temperature,args.thermostat,params[0], np.sqrt(cov[0,0]), params[1], np.sqrt(cov[1,1]) ))
out_short.close()

#
# Long-time correlation function
#
base,extension=path.splitext(path.basename(args.corrname))
np.savetxt(base+'_long'+extension, 
			np.column_stack((times, KmctLong, KmctLong/C[0])), 
			fmt=['%.14g','%.14g','%.14g'], 
			header='time Klong Klong/K[0]')


