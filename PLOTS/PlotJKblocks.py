#!/usr/bin/env python

import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt



parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
# parser.add_argument('-f','--filename', default="", help='name of the file you want to plot')
parser.add_argument('-T','--T', default='5.0', help='Temperature')
parser.add_argument('-N','--N', default=1080, help='number of particles')
parser.add_argument('-t','--thermostat', default='NVT', help='thermostat')
parser.add_argument('-o','--obs', choices=['CPP','CFP','CFF','K'], default='CPP', help='What to plot')
parser.add_argument('-M','--M', default=3, help='M of Gaver-Stehfest (only for K)')
parser.add_argument('--nomean', action='store_true', help='if activated, mean is not plotted')
args = parser.parse_args()


# LOAD DATA JK
path='../OUTPUT/T{}/N{}/'.format(args.T, args.N)

times=np.load(path+'times_{}.npy'.format(args.thermostat))
if args.obs=='CPP' or args.obs=='CFP' or args.obs=='CFF':
	filename=path+'/{}JK_{}.npy'.format(args.obs, args.thermostat)
	corr=np.load(filename).item()
elif args.obs=='K':
	filename=path+'/noisecorrJK_{}_M{}.npy'.format(args.thermostat,args.M)
	corr=np.load(filename).item()['combine']
	av=np.loadtxt(path+'noisecorr_{}_combine_M{}.txt'.format(args.thermostat,args.M))

# LOAD DATA AVERAGE (average is different than mean because of the order of the operations)


#PLOT
fig=plt.subplot(111, xscale='log')
if args.obs=='K': 
	plt.ylim(top=1.1*corr['mean'][0])
plt.title(args.obs+', T='+args.T+'  (all JK blocks + mean)')
nblo=len(corr['blocksJK'])
for iblo in range(nblo):
	plt.plot(times,corr['blocksJK'][iblo])
if not args.nomean:
	plt.errorbar(times, corr['mean'], yerr=corr['errJK'], linewidth=3, color='black')
try:
	plt.plot(av[:,0],av[:,1], label="Average", color='grey', linewidth=3)
except NameError:
	pass

plt.show()