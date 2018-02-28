#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt

nsamples=4
T=0.466
data={'samples':[],'mean':[],'err':[]}
data['samples']=[np.loadtxt('../OUTPUT/T'+str(T)+'/N1080/S'+str(i)+'/corrF.txt',comments='#',usecols=(1,4,5), dtype=np.float64) for i in range(nsamples)]
data['mean']=np.mean(data['samples'],axis=0)
data['err']=np.std(data['samples'],axis=0)/np.sqrt(nsamples-1)


plt.xlabel('$t$')
plt.ylabel('$C(t)/C(0)$')
plt.xscale('log')

for i in range(nsamples):
	plt.plot(data['samples'][i][:,0], data['samples'][i][:,1], linewidth='1.0', color='magenta')
	plt.plot(data['samples'][i][:,0], data['samples'][i][:,2], linewidth='1.0', color='orange')

plt.errorbar(data['mean'][:,0], data['mean'][:,1], data['err'][:,1], color='blue')
plt.errorbar(data['mean'][:,0], data['mean'][:,2], data['err'][:,2], color='red')
plt.grid()
plt.legend(loc='upper right')
plt.show()
