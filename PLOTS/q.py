#!/usr/bin/env python

import numpy as np
import sys
from matplotlib import pyplot as plt
from glob import glob
from scipy.stats import sem

temperatures=[5.0, 2.0, 1.0, 0.8, 0.7, 0.6, 0.55, 0.52, 0.49]
# temperatures=[5.0]
nT=len(temperatures)


ovlps={}
for iT in range(nT):

	directorio='/storage4Tb/STRUCTURAL-GLASS/OUTPUT/T{}/N65/shift/'.format(temperatures[iT])
	filenames=glob(directorio+'/S*/OverlapTrajectoryISirep*.txt')

	t=np.loadtxt(filenames[0], usecols=(1))
	allq,allqIS = np.ndarray( len(t) ), np.ndarray( len(t) )

	array=np.array([np.loadtxt(f, usecols=(4,5)) for f in filenames])
	q  =array[:,:,0]
	qIS=array[:,:,1]

	ovlps  [temperatures[iT]] = {	
									'time':t, 
									'term': {'mean': q.mean(  axis=0), 'err':sem(q,  axis=0)}, 
									'IS': {'mean': qIS.mean(axis=0), 'err':sem(qIS,axis=0)}
								}







dt=0.0025
fig = plt.figure()
colors=['#ed0000','#d70216','#c1042c','#ac0642','#ac0642','#ac0642','#6b0d85','#560f9b','#4011b2','#2b13c8','#1515de','#0018f5']

for iT in range(nT):
	ax1 = plt.subplot(3, 3, iT+1)
	ax1.errorbar(ovlps[temperatures[iT]]['time']*dt, ovlps[temperatures[iT]]['term']['mean'], yerr=ovlps[temperatures[iT]]['term']['err'], linestyle=':', color='red')
	ax1.errorbar(ovlps[temperatures[iT]]['time']*dt, ovlps[temperatures[iT]][ 'IS' ]['mean'], yerr=ovlps[temperatures[iT]][ 'IS' ]['err'], label=r'$T={}$'.format(temperatures[iT]), linestyle='-', color='black')
	ax1.set_xscale('log')
	ax1.set_xlabel(r'$t$')
	ax1.set_ylabel(r'$q$')
	ax1.legend()
plt.savefig('./FIGURES/OverlapTrajectoryIS.png')
plt.show()
sys.exit()

