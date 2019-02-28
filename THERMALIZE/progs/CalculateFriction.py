#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import lib.module_measurements as med

# READ COMMAND-LINE ARGUMENTS
parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
parser.add_argument('corrname', help='correlation function we want to read')
parser.add_argument('--temperature', type=float, required=True, help='temperature is needed to transform the correlator into friction')
parser.add_argument('--density', type=float, required=False, default=1.2, help='density is optional because I always work at rho=1.2')
parser.add_argument('--tstar', type=float, required=False, default=None, help='no integration for t>tstar')
args = parser.parse_args()

if (not args.tstar==None) and args.tstar<=0: raise ValueError('tstar should be positive ({})'.format(args.tstar))

data=np.loadtxt(args.corrname)
times=data[:,0]
C=data[:,1]
# errC=data[:,2]
del data


istar= -1 if args.tstar==None else np.where(times>args.tstar)[0][0]


friction = np.trapz(C[:istar], x=times[:istar]) * args.density/args.temperature
print(friction)


# # DRAW THE CORRELATION
# from matplotlib import pyplot as plt
# plt.subplot(111, xscale="log", yscale="linear")
# # # plt.errorbar(times, C, yerr=errC)
# plt.plot(times, C)
# plt.grid()
# plt.show()

# def f(x, a1, a2):
# 	return a1/np.cosh(a2*x)


# ishort=np.where(times>0.05)[0][0]

# params=[C[0],1]
# params,dummy=curve_fit(f, times[:ishort], C[:ishort], p0=params)
# plt.subplot(111, xscale="log", yscale="linear", ylim=(0,1+C[0]))
# plt.plot(times, C)
# plt.plot(times, f(times,params[0],params[1]))
# plt.grid()
# plt.show()


