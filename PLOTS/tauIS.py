#!/usr/bin/env python


import numpy as np
from matplotlib import pyplot as plt
from lib import module_histo as his
from scipy.stats import sem
from scipy.optimize import curve_fit
import pandas as pd
import sys
from glob import glob

dt=0.0025
deltaE=1e-5

temperatures=[0.49, 0.52, 0.55, 0.6, 0.7, 0.8, 1.0]
# temperatures=[0.6]
nT=len(temperatures)

# /storage4Tb/STRUCTURAL-GLASS/OUTPUT/T0.6/N65/shift/S0/chunksIS/


tauIS={}
for iT in range(nT):
	print('T = ',temperatures[iT])
	directorio='/storage4Tb/STRUCTURAL-GLASS/OUTPUT/T{}/N65/shift/'.format(temperatures[iT])

	# 
	# INHERENT STRUCTURES
	# 
	filenames=glob(directorio+'S*/chunksIS/MinimizeSegment.txt')
	arrays=[np.loadtxt(f)[:,2] for f in filenames]

	taus=[]
	for iarray in range(len(arrays)):
		tempo=0
		for i in range(1,len(arrays[iarray])):
			# print(i, arrays[iarray][i])
			
			if np.isclose(arrays[iarray][i],arrays[iarray][i-1] , atol=1e-5):
				tempo+=1
			else:
				taus.append(tempo)
				tempo=0

	taus=np.array(taus)+1e-6 #We add a small number, to push the distribution away from zero and plot it in logs
	nbin=max(int(len(taus)/40), 10)
	histotau = his.logbinning(taus, n=nbin, xmin=np.min(taus), xmax=-1, silent=True) #
	tauIS[temperatures[iT]] = {'taus': np.array(taus), 'histotau':histotau, 'mean': np.mean(taus), 'err': sem(taus)}


# 
# Plot distributions
# 
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\tau_\mathrm{IS}$')
ax1.set_ylabel(r'$h(\tau_\mathrm{IS})$')

for iT in range(len(temperatures)):
	T=temperatures[iT]
	ax1.plot(tauIS[T]['histotau'][:,1]*dt, tauIS[T]['histotau'][:,2], label=r'$T = {}$'.format(T) )

plt.show()


# 
# Plot mean values
# 
tis=np.array([tauIS[temperatures[iT]]['mean'] for iT in range(nT)])*dt
tiserr=np.array([tauIS[temperatures[iT]]['err'] for iT in range(nT)])*dt
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_xscale('linear')
ax1.set_xlabel(r'$T$')
ax1.set_ylabel(r'$\tau_\mathrm{IS}$')
ax1.errorbar(temperatures, tis, yerr=tiserr, label=r'$T = {}$'.format(T) )
plt.show()


np.savetxt('../THERMALIZE/data/tauIS.txt' ,np.column_stack( [temperatures, tis, tiserr]), fmt='%g %.10f %g' )
