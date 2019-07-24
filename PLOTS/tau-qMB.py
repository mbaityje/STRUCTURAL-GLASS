#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from lib import module_histo as his
from scipy.stats import sem
from scipy.optimize import curve_fit
import pandas as pd
import sys

dt=0.0025
ta = pd.read_csv('../THERMALIZE/data/tau_alpha.txt',sep=' ', comment='#')

def f(x, A, Tc, eta):
	return A/(x-Tc)**eta
p0=(1,0.43,1.5)
imin=6#4
imax=len(ta['tau'])+1
popt, pcov = curve_fit(f, ta['T'].values[imin:imax], ta['tau'].values[imin:imax], p0=p0, bounds=((0,-np.inf,0),(np.inf,0.459,5)) )
plt.xlabel('T-T^*')
plt.ylabel(r'$\tau_\alpha$')
plt.yscale('log')
plt.xscale('log')
plt.text(0.1,10,'$T^*=${}'.format(popt[1]))
plt.plot(ta['T']-popt[1],ta['tau'], 'o-')
plt.plot(ta['T']-popt[1],f(ta['T'].values, popt[0], popt[1], popt[2]))
# plt.savefig('taualphaTN65.png')
plt.show()

# Inherent structure time
tIS=np.loadtxt('../THERMALIZE/data/tauIS.txt')


temperatures=[0.49, 0.52, 0.55, 0.6, 0.7, 0.8, 1.0]

obs={}
for T in temperatures:

	directorio='/storage4Tb/STRUCTURAL-GLASS/OUTPUT/T{}/N65/shift/'.format(T)

	#
	# ANALYZE TAU
	#


	tauFILE=directorio+'/tauMB.txt'
	# print(tauFILE)

	taus=np.loadtxt(tauFILE)*dt #+1e-6 #Add 1e-6 to avoid having zeros
	# nbin=max(int(len(taus)/100), 10)
	nbin=30
	taumedian=np.median(taus)
	taumean  = taus.mean()
	histotau = his.logbinning(taus, n=nbin, xmin=dt, xmax=-1, silent=True) #xmin=taus.min()
	taualpha = ta.loc[ (ta['T']==T) ]['tau'].values.item()
	print('T={}\ttau= {}\taveragetauMB= {}\ttaumedianMB= {}'.format(T, taualpha, taumean, taumedian*dt ))

	tauIS=tIS[np.where(tIS[:,0]==T)[0][0]][1]
	print('tauIS:',tauIS)

	#
	# ANALYZE OVERLAP
	#

	samples=np.arange(10)
	ovlp=[]
	for isam in samples:
		qFILE=directorio+'/S{}/chunksIS/qMB.txt'.format(isam)
		try:
			ovlp.append(np.loadtxt(qFILE)) # q(n), qintra, qextra
		except:
			pass

	qintra    = np.mean([ovlp[isam][0][1] for isam in range(len(ovlp))])
	qintraErr =     sem([ovlp[isam][0][1] for isam in range(len(ovlp))])
	qextra    = np.mean([ovlp[isam][0][2] for isam in range(len(ovlp))])
	qextraErr =     sem([ovlp[isam][0][2] for isam in range(len(ovlp))])
	print('T={}\tqintra= {}\tqextra= {}'.format(T, qintra, qextra ))


	# maxdist = np.min([len(ovlp[i]) for i in range(len(ovlp))])
	# qMBevol    = np.array([np.mean([ovlp[isam][dist][0] for isam in range(len(ovlp))]) for dist in range(maxdist)])
	# qMBevolErr = np.array([    sem([ovlp[isam][dist][0] for isam in range(len(ovlp))]) for dist in range(maxdist)])
	maxdist = np.max([len(ovlp[i]) for i in range(len(ovlp))])
	qMBevol    = np.array([  np.nanmean( [(ovlp[isam][dist][0] if len(ovlp[isam])>dist else np.nan) for isam in range(len(ovlp))]  ) for dist in range(maxdist)])
	qMBevolErr = np.array([  sem( [(ovlp[isam][dist][0] if len(ovlp[isam])>dist else np.nan) for isam in range(len(ovlp))], nan_policy='omit') for dist in range(maxdist)])


	#Record observables of this temperature
	obs[T]= {
					'histotau' : {'mean': histotau}, \
					'taualpha' : {'mean': taualpha}, \
					'taumean'  : {'mean': taumean,  'err':sem(taus)}, \
					'taumedian': {'mean': taumedian}, \
					'qintra'   : {'mean': qintra,  'err':qintraErr}, \
					'qextra'   : {'mean': qextra,  'err':qextraErr}, \
					'qMBevol'  : {'mean': qMBevol, 'err':qMBevolErr} \
					}


	#Plot 
	fig = plt.figure()
	fig.suptitle('T = {}'.format(T), fontsize=11)
	ax1 = fig.add_subplot(1, 1, 1)
	# ax1 = fig.add_subplot(2, 1, 1)
	ax1.set_yscale('log')
	ax1.set_xscale('log')
	ax1.set_xlabel(r'$\tau_\mathrm{MB}$')
	ax1.set_ylabel(r'$P(\tau_\mathrm{MB})$')
	plt.axvline(x=taualpha, color='red', label=r'$\langle\tau\rangle$')
	plt.axvline(x=taumean, color='grey', label=r'$\langle\tau_\mathrm{MB}\rangle$')
	plt.axvline(x=taumedian, color='black', label=r'$\tau_\mathrm{MB,1/2}$')
	# plt.plot(histotau[:,1], histotau[:,2], label=r'$h(\tau_\mathrm{MB})$')
	plt.plot(histotau[:,1], histotau[:,3], marker='x', label=r'$h(\tau_\mathrm{MB})$')
	plt.legend(loc='best')
	# ax2 = fig.add_subplot(2, 1, 2)
	# ax2.set_ylabel(r'$q_\mathrm{MB}$')
	# ax2.set_xlabel('basin distance')
	# ax2.set_ylim(top=1)
	# plt.plot(range(maxdist), np.full(maxdist,qintra), color='grey', label=r'$q_\mathrm{intra}$')
	# plt.plot(range(maxdist), np.full(maxdist,qextra), color='black', label=r'$q_\mathrm{extra}$')
	# plt.errorbar(range(maxdist), qMBevol, yerr=qMBevolErr, label=r'$q_\mathrm{MB}$')
	# plt.legend(loc='best')
	# plt.savefig('./FIGURES/tau-qMB_T{}.png'.format(T))

	# extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	# fig.savefig('./FIGURES/histotau_T{}.png'.format(T), bbox_inches=extent.expanded(1.5, 1.5))
	fig.savefig('./FIGURES/histotau_T{}.png'.format(T))
	plt.show()



# Plots comparing different temperatures
fig = plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)

# ax4 = fig.add_subplot(4, 1, 4)

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$\tau_\mathrm{MB}$')
ax1.set_ylabel(r'$P(\tau_\mathrm{MB})$')
ax1.set_xlim(left=dt, right=1e5)
ax2.set_ylabel(r'$q_\mathrm{MB}$')
ax2.set_xlabel(r'basin distance')
ax2.set_ylim(top=1)
ax2.set_xlim(right=500)

# Tc=popt[1]
# eta=popt[2]
# ax4.set_ylabel(r'$q_\mathrm{MB}$')
# ax4.set_xlabel(r'$T$')

ax3.set_ylabel(r'$q_\mathrm{intra}$')
ax3.set_xlabel(r'$T$')

for T in temperatures:
	# ax1.plot(obs[T]['histotau']['mean'][:,1], obs[T]['histotau']['mean'][:,2], marker='x', label=r'$T = {}$'.format(T) )
	ax1.plot(obs[T]['histotau']['mean'][:,1], obs[T]['histotau']['mean'][:,3], marker='x', label=r'$T = {}$'.format(T) )
	ax2.errorbar(range(len(obs[T]['qMBevol']['mean'])), obs[T]['qMBevol']['mean'], yerr=obs[T]['qMBevol']['err'], label=r'$T = {}$'.format(T))

	# ax4.errorbar((np.arange(len(obs[T]['qMBevol']['mean']))-Tc)**(-eta), obs[T]['qMBevol']['mean'], yerr=obs[T]['qMBevol']['err'], label=r'$T = {}$'.format(T))

ax3.errorbar(temperatures, [obs[T]['qintra']['mean'] for T in temperatures], yerr=[obs[T]['qintra']['err'] for T in temperatures])
ax2.legend()
plt.show()





#
# PLOTS of tauIS/tauMB vs T
# 
tauMB=np.array([obs[T]['taumean']['mean'] for T in temperatures])
tauMBerr=np.array([obs[T]['taumean']['err'] for T in temperatures])
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_yscale('linear')
ax1.set_xscale('linear')
ax1.set_xlabel(r'$T$')
ax1.set_ylabel(r'$\frac{\langle\tau_\mathrm{IS}\rangle}{\langle\tau_\mathrm{MB}\rangle}$')
ax1.errorbar(temperatures, tIS[:,1]/tauMB, yerr=(tIS[:,1]*tauMBerr+tIS[:,2]*tauMB)/np.square(tauMB), marker='o', color='darkred')
# ax1.errorbar(temperatures, tIS[:,1]/tauMB, yerr=(tIS[:,2])/tauMB )
# ax1.errorbar(temperatures, tIS[:,1]/tauMB, yerr=(tIS[:,1]*tauMBerr)/(tauMB*tauMB) )
plt.savefig('./FIGURES/tauIS-on-tauMB.png'.format(T))
plt.show()
