#!/usr/bin/env python

import sys
import numpy as np
import argparse
from scipy.interpolate import interp1d
from scipy import integrate
import scipy
import lib.module_measurements as med
from matplotlib import pyplot as plt
from lib.beylkin import Beylkin


def FitFunction(x, y, ncoef=10):
	'''Fit function as a sum of exponential functions'''
	b = Beylkin(decaying=True)
	b.driver_load(x, y, ncoef)
	return b.prony_function


	# checks.IntegrateCFP()
	# checks.PlotCPPcheck()
	# checks.CalculateCFPregularized()
	# checks.PlotCFPcheck()


class CorrelationConsistency:
	def __init__(self):
		self.ReadArgs()
		self.PrintArgs()
		self.SetShowPlots(self.args.showplots)
		self.ReadCorrelations()
		self.CalculateCPPdotJK()
		if self.showplots: 
			self.PlotK()
			self.PlotCorrelations()
			self.PlotCPP()
		self.CalculateCFPJK()
		if self.showplots: self.PlotCFPcheck()
		self.IntegrateCFPJK()
		self.PlotCPPcheck()
		return

	def SetShowPlots(self, switch):
		assert(isinstance(switch, bool))
		self.showplots = switch

	def ReadArgs(self):	
		''' Read command-line arguments '''
		parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
		parser.add_argument('filenameK', help='Filename with the JK blocks dictionary in binary (e.g. noisecorrJK_NVT_M4.npy)')
		parser.add_argument('--thermostat', required=False, default='NVE', choices=['NVE','NVT'], help='thermostat, necessary for filenames')
		parser.add_argument('--showplots', action='store_true', help='If flag is activated shows plots of the functions')
		self.args = parser.parse_args()
		print(sys.argv)

	def PrintArgs(self):
		print('filenameK:',self.args.filenameK)
		print('thermostat:',self.args.thermostat)


	def ReadCorrelations(self):
		self.CFF  = np.load('CFFJK_'+str(self.args.thermostat)+'.npy') # access as CFF.item()['mean']
		self.CFP  = np.load('CFPJK_'+str(self.args.thermostat)+'.npy')
		self.CPP  = np.load('CPPJK_'+str(self.args.thermostat)+'.npy')
		self.K    = np.load(self.args.filenameK)
		self.nblo = len(self.K.item()['combine']['blocksJK'])
		self.times= np.load('times_'+str(self.args.thermostat)+'.npy')
		self.nt   = len(self.times)


		if not (len(self.CFP.item()['mean'])==self.nt and len(self.CFF.item()['mean'])==self.nt and len(self.CPP.item()['mean'])==self.nt):
			raise IndexError('The number of correlations ('+str(len(self.CFP.item()['mean']))+' or '+str(len(self.CFF.item()['mean']))+' or '+str(len(self.CPP.item()['mean']))+') is inconsistent with the number of times('+str(self.nt)+')')
		if len(self.K.item()['combine']['mean']) != self.nt:
			raise ValueError('Number of times in K is different from the one in self.times')

	def PlotK(self):
		fig, ax = plt.subplots(1,1)
		ax.set_xscale('log')
		ax.set_xlabel('$t$')
		ax.set_ylabel('$\mathcal{K}(t)$')
		for iblo in range(self.nblo):
			ax.plot(self.times, self.K.item()['combine']['blocksJK'][iblo], '.', label='$\mathcal{K}$')
		ax.plot(self.times, self.K.item()['combine']['mean'], label='$\mathcal{K}$', color='black')
		ax.legend()
		plt.show()

	def PlotCorrelations(self):
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1)
		ax1.set_ylabel('$C^{FP}(t)$')
		ax1.get_xaxis().set_visible(False)
		ax1.set_xscale('log')
		ax1.plot(self.times, self.CFP.item()['mean'])
		ax2.set_ylabel('$C^{F}(t)$')
		ax2.get_xaxis().set_visible(False)
		ax2.set_xscale('log')
		ax2.plot(self.times, self.CFF.item()['mean'])
		ax3.set_ylabel('$C^{P}(t)$')
		ax3.get_xaxis().set_visible(False)
		ax3.set_xscale('log')
		ax3.plot(self.times, self.CPP.item()['mean'])
		ax4.set_ylabel('$\mathcal{K}(t)$')
		ax4.set_xlabel('$t$')
		ax4.get_xaxis().set_visible(True)
		ax4.set_xscale('log')
		ax4.plot(self.times, self.K.item()['combine']['mean'], color='darkred')
		plt.show()

	def PlotCPP(self):
		fig, ax = plt.subplots(1,1)
		ax.set_ylabel('Correlation Function')
		ax.set_xlabel('$t$')
		ax.set_xscale('log')
		ax.plot(self.times, self.CPP.item()['mean'], color='darkblue' , label='$C^{PP}(t)$')
		ax.errorbar(self.times, self.CPPdot['mean'], yerr=self.CPPdot['errJK'], color='darkgreen', label='$\dot{C}^{PP}(t)$')
		ax.plot(self.times, -self.CFP.item()['mean'], color='darkred'  , label='$-C^{FP}(t) = C^{PF}(t)$')
		ax.legend()
		plt.show()


	def CalculateCPPdotJK(self):
		self.CPPdot={
			'mean':np.ndarray(self.nt, dtype=np.float64),
			'blocksJK':np.ndarray((self.nblo, self.nt), dtype=np.float64),
			'errJK':np.ndarray(self.nt, dtype=np.float64),
			}
		for iblo in range(self.nblo):
			self.CPPdot['blocksJK'][iblo]=self.CalculateDerivative(self.CPP.item()['blocksJK'][iblo])
		self.CPPdot['mean' ] = self.CPPdot['blocksJK'].mean(axis=0)
		self.CPPdot['errJK'] = np.sqrt((self.nblo-1)*(np.square(self.CPPdot['blocksJK']).mean(axis=0) - np.square(self.CPPdot['mean']) ) )



	def CalculateDerivative(self, f):
		fdot=np.ndarray(self.nt, dtype=np.float64)
		t=self.times
		dfp=np.float64(f[1]-f[0])
		dtp=np.float64(t[1]-t[0])
		fdot[0]=dfp/dtp
		for i in range(1,self.nt-1):
			dfm=dfp
			dtm=dtp
			dfp=f[i+1]-f[i]
			dtp=t[i+1]-t[i]
			fdot[i]=0.5*(dfp/dtp+dfm/dtm)
		fdot[self.nt-1] = dfm/dtm
		return fdot

	def CalculateCFPJK(self):
		self.CFPcheck={
			'mean':np.ndarray(self.nt, dtype=np.float64),
			'blocksJK':np.ndarray((self.nblo, self.nt), dtype=np.float64),
			'errJK':np.ndarray(self.nt, dtype=np.float64),
			}
		for iblo in range(self.nblo):
			self.CFPcheck['blocksJK'][iblo]=self.CalculateCFP(self.CPP.item()['blocksJK'][iblo], self.K.item()['combine']['blocksJK'][iblo])
		self.CFPcheck['mean' ] = self.CFPcheck['blocksJK'].mean(axis=0)
		self.CFPcheck['errJK'] = np.sqrt((self.nblo-1)*(np.square(self.CFPcheck['blocksJK']).mean(axis=0) - np.square(self.CFPcheck['mean']) ) )

	def CalculateCFP(self, cpp, K):
		'''
		The noise correlation K is the memory kernel of CPF=-CFP.
		This function calculates CFP using K, in order to compare it to the measured CFP.
		'''
		t=self.times
		interpK=interp1d(t, K, kind='cubic')

		def trapeze(i):
			temp=0
			for j in range(1, i+1):
				# print('i:%d  j:%d  ti-tj: %f  ti-tjm1: %f'%(i,j,t[i]-t[j],t[i]-t[j-1]))
				temp += ( interpK(t[i]-t[j])*cpp[j] + interpK(t[i]-t[j-1])*cpp[j-1] )*(t[j]-t[j-1])
			return 0.5*temp 

		tempCFP=np.ndarray(self.nt ,dtype=np.float64)
		tempCFP[0]=0
		for i in range(0,self.nt):
			print('\rtrapeze iteration', i, end='')
			tempCFP[i]=trapeze(i) #in this version of the program we assume that the input data already was divided by kT
		print('')
		return tempCFP

	def IntegrateCFPJK(self):
		self.CPPcheck={
			'mean':np.ndarray(self.nt, dtype=np.float64),
			'blocksJK':np.ndarray((self.nblo, self.nt), dtype=np.float64),
			'errJK':np.ndarray(self.nt, dtype=np.float64),
			}
		for iblo in range(self.nblo):
			self.CPPcheck['blocksJK'][iblo]=self.IntegrateCFP(self.CPP.item()['blocksJK'][iblo][0], self.CFPcheck['blocksJK'][iblo])
		self.CPPcheck['mean' ] = self.CPPcheck['blocksJK'].mean(axis=0)
		self.CPPcheck['errJK'] = np.sqrt((self.nblo-1)*(np.square(self.CPPcheck['blocksJK']).mean(axis=0) - np.square(self.CPPcheck['mean']) ) )


	def IntegrateCFP(self, cpp0, cfpcheck):
		CPPcheck=np.full(self.nt, cpp0)
		for i in range(self.nt):
			CPPcheck[i]-=integrate.trapz(cfpcheck[:i+1], x=self.times[:i+1])
		return CPPcheck

	def PlotCFPcheck(self):
		fig,ax=plt.subplots(1,1)
		ax.set_xscale('log')
		ax.errorbar(self.times, self.CFPcheck['mean'], yerr=self.CFPcheck['errJK'], color='red'    , label='$C^{FP}$ From memory kernel')
		ax.errorbar(self.times, self.CFP.item()['mean'], yerr=self.CFP.item()['errJK'], color='darkred', label='$C^{FP}$ Direct measurement')
		ax.errorbar(self.times, self.K.item()['combine']['mean']*np.max(self.CFPcheck['mean'])/self.K.item()['combine']['mean'][0], yerr=self.K.item()['combine']['errJK']*np.max(self.CFPcheck['mean'])/self.K.item()['combine']['mean'][0], color='grey', label='$\mathcal{K}$')
		ax.errorbar(self.times, self.CPP.item()['mean']*np.max(self.CFPcheck['mean'])/self.CPP.item()['mean'][0], yerr=self.CPP.item()['errJK']*np.max(self.CFPcheck['mean'])/self.CPP.item()['mean'][0], color='lightgrey', label='$C^{PP}$')

		ax.set_ylabel('$C(t)$')
		ax.set_xlabel('$t$')
		ax.legend()
		plt.show()

	def PlotCPPcheck(self):
		fig,ax=plt.subplots(1,1)
		ax.set_xscale('log')
		ax.errorbar(self.times, self.CPPcheck['mean'  ], yerr=self.CPPcheck['errJK'  ], color='red'    , label='$C^{PP}(t)$ From memory kernel')
		ax.errorbar(self.times, self.CPP.item()['mean'], yerr=self.CPP.item()['errJK'], color='darkred', label='$C^{PP}(t)$ Direct measurement')
		ax.set_ylabel('$C(t)$')
		ax.set_xlabel('$t$')
		ax.set_ylim((-1,self.CPP.item()['mean'].max()))
		ax.legend()
		plt.show()


	def WriteCPPcheck(self):
		np.savetxt('CPPcheck_{}_M{}.txt'.format(args.thermostat,args.M), 
			np.column_stack((
				self.times, 
				self.CPPcheck['mean'  ], self.CPPcheck['errJK'  ],
				self.CPP.item()['mean'], self.CPP.item()['errJK']
				)),
			fmt=['%g','%.14g','%.14g','%.14g','%.14g'],
			header='time CPPcheck errCPPcheck CPP errCPP'
			)

if __name__ == "__main__":
	checks=CorrelationConsistency()
	checks.WriteCPPcheck()










