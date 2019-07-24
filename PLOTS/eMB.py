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
temperatures=[0.49, 0.52, 0.55, 0.6, 0.7, 0.8, 1.0]
nT=len(temperatures)

E={ 'E':        np.ndarray((nT,2)),
	'Eis':      np.ndarray((nT,2)),
	'Eridge':   np.ndarray((nT,2)),
	'Emb':      np.ndarray((nT,2)),
	'EridgeMB': np.ndarray((nT,2))}

for iT in range(nT):

	directorio='/storage4Tb/STRUCTURAL-GLASS/OUTPUT/T{}/N65/shift/'.format(temperatures[iT])

	# 
	# INHERENT STRUCTURES
	# 
	filenames=glob(directorio+'S*/chunksIS/elistIS.txt')
	arrays=[np.loadtxt(f)[:,1] for f in filenames]
	enersIS=np.concatenate(arrays)

	filenames=glob(directorio+'S*/chunksIS/elistRidge.txt')
	arrays=[np.loadtxt(f)[:,1] for f in filenames]
	enersRidge=np.concatenate(arrays)

	E['Eis'   ][iT] = [enersIS.mean(), sem(enersIS)]
	E['Eridge'][iT] = [enersRidge.mean(), sem(enersRidge)]


	#
	# METABASINS
	#

	embFILE=directorio+'/EMB.txt'
	enerRidgeFILE=directorio+'/EridgeMB.txt'

	enersMB=np.loadtxt(embFILE)
	enersRidgeMB=np.loadtxt(enerRidgeFILE)[:,1]

	E['Emb'     ][iT] = [enersMB.mean(), sem(enersMB)]
	E['EridgeMB'][iT] = [enersRidgeMB.mean(), sem(enersRidgeMB)]

	# print('T={}\tenerMB= {}±{}\tenerRidgeMB= {}±{}'.format(temperatures[iT], E['Emb'][iT][0], E['Emb'][iT][1], E['EridgeMB'][iT][0], E['EridgeMB'][iT][1]  ) )



	#
	# AVERAGE ENERGY
	# 
	filenames=glob(directorio+'S*/_aftergap_shift_NVT.txt')
	eners=[np.loadtxt(f, comments='#')[:,3].mean() for f in filenames]

	E['E'][iT] = [np.mean(eners), sem(eners)]

#
# Plot comparing different temperatures
#
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

ax1.set_yscale('linear')
ax1.set_xscale('linear')
ax1.set_xlim((0,1))
ax1.set_xlabel(r'$T$')
ax1.set_ylabel(r'$E$')

ax1.errorbar(temperatures, E['E'][:,0], yerr=E['E'][:,1], label=r'$E$', linestyle=':', color='black')

ax1.errorbar(temperatures, E['Eis'][:,0], yerr=E['Eis'][:,1], label=r'$E_\mathrm{IS}$', linestyle='-', color='red')
ax1.errorbar(temperatures, E['Eridge'][:,0], yerr=E['Eridge'][:,1], label=r'$E^\mathrm{ridge}_\mathrm{IS}$', linestyle='--', color='red')

ax1.errorbar(temperatures, E['Emb'][:,0], yerr=E['Emb'][:,1], label=r'$E_\mathrm{MB}$', linestyle='-', color='darkgreen')
ax1.errorbar(temperatures, E['EridgeMB'][:,0], yerr=E['EridgeMB'][:,1], label=r'$E^\mathrm{ridge}_\mathrm{MB}$' , linestyle='--', color='darkgreen')
ax1.legend()

# Inset
left, bottom, width, height = [0.2, 0.2, 0.3, 0.3]
axi = fig.add_axes([left, bottom, width, height])
axi.errorbar(temperatures, E['Eis'][:,0], yerr=E['Eis'][:,1], label=r'$E_\mathrm{IS}$', linestyle='-', color='red')
axi.errorbar(temperatures, E['Eridge'][:,0], yerr=E['Eridge'][:,1], label=r'$E^\mathrm{ridge}_\mathrm{IS}$', linestyle='--', color='red')
axi.errorbar(temperatures, E['Emb'][:,0], yerr=E['Emb'][:,1], label=r'$E_\mathrm{MB}$', linestyle='-', color='darkgreen')
axi.errorbar(temperatures, E['EridgeMB'][:,0], yerr=E['EridgeMB'][:,1], label=r'$E^\mathrm{ridge}_\mathrm{MB}$' , linestyle='--', color='darkgreen')
plt.savefig('./FIGURES/eMB.png')
plt.show()




