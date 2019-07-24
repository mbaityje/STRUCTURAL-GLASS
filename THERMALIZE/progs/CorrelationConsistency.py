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

class CorrelationConsistency:
	def __init__(self):
		return

	def ReadArgs(self):	
		''' Read command-line arguments '''
		parser = argparse.ArgumentParser(prog='python '+sys.argv[0]+' [--hoomd-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"', add_help=True)
		parser.add_argument('-N','--Natoms', type=int, required=True, help='Number of particles')
		parser.add_argument('-M','--M', type=int, required=False, default=5, help='(Half) Number of Gaver-Stehfest coefficients')
		parser.add_argument('-L','--L', type=np.float64, required=True, help='Box Size (assumed cubic)')
		parser.add_argument('-T','--temperature', type=float, required=True, help='Temperature')
		parser.add_argument('--tstar', type=float, required=False, default=1.0, help='Maximum time for integral of CPP(t)')
		parser.add_argument('--thermostat', required=False, default='NVE', choices=['NVE','NVT'], help='thermostat, necessary for filenames')
		parser.add_argument('--smoothK', required=False, default=None, choices=['None','ReLU','linear','exp','custom'], help='Apply smoothing to the K(t)')
		self.args = parser.parse_args()
		if self.args.smoothK=='None': self.args.smoothK=None
		self.CheckArgs()
		print(sys.argv)

	def PrintArgs(self):
		print('N:',self.args.Natoms)
		print('L:',self.args.L)
		print('T:',self.args.temperature)
		print('thermostat:',self.args.thermostat)

	def CheckArgs(self):
		if not self.args.Natoms>0: raise ValueError('L ({}) must be positive'.format(args.Natoms))
		if not self.args.L>0: raise ValueError('L ({}) must be positive'.format(args.L))
		if self.args.temperature<0:
			raise ValueError('Temperature='+str(self.args.temperature)+'. Must be positive. Aborting.')
		else:
			self.invT=np.float64(1./self.args.temperature)

	def ReadCorrelations(self):
		self.CFF=np.load('CFF_'+str(self.args.thermostat)+'.npy') # access as CFP.item()['mean']
		self.CFP=np.load('CFP_'+str(self.args.thermostat)+'.npy')
		self.CPP=np.load('CPPJK_'+str(self.args.thermostat)+'.npy')
		self.times=np.load('times_'+str(self.args.thermostat)+'.npy')
		self.nt=len(self.times)
		readK=np.loadtxt('noisecorr_'+self.args.thermostat+'_combine_M'+str(self.args.M)+'.txt')[:,0:2]
		self.Kread=readK[:,1]*self.args.temperature #at some point I changed the def of K by a factor 1/T, so sometimes we need to remultiply it back
		self.SmoothK(action=self.args.smoothK)


		if not (len(self.CFP.item()['mean'])==self.nt and len(self.CFF.item()['mean'])==self.nt and len(self.CPP.item()['mean'])==self.nt):
			raise IndexError('The number of correlations ('+str(len(self.CFP.item()['mean']))+' or '+str(len(self.CFF.item()['mean']))+' or '+str(len(self.CPP.item()['mean']))+') is inconsistent with the number of times('+str(self.nt)+')')
		if len(readK[:,0]) != self.nt:
			raise ValueError('Number of times in K is different from the one in self.times')
		if any(np.abs(readK[:,0]-self.times)>1e-8):
			raise ValueError('Times in K and in self.times are inconsistent')

	def IdentifyMaxK(self):
		#First find the min
		for i in range(1,self.nt):
			if self.Kread[i]>self.Kread[i-1]:
				imin=i-1
				break
		#Then find the local maximum
		for i in range(imin+1,self.nt):
			if self.Kread[i]<self.Kread[i-1]:
				imax=i-1
				break
		#Then find the dip
		for i in range(imax+1,self.nt):
			if self.Kread[i]>self.Kread[i-1] and self.Kread[i]<0:
				idip=i-1
				break
		#The find when it cuts zero
		for i in range(idip+1,self.nt):
			if self.Kread[i]>0:
				izero=i
				break
		print('times[imin] =', self.times[imin],'\ttimes[imax] =', self.times[imax], '\ttimes[idip] =', self.times[idip],'\ttimes[izero] =', self.times[izero])
		return imin,imax,idip,izero

	def SmoothK(self, action=None):

		print('action = ',action)
		if action==None:
			self.K=self.Kread.copy()
			return

		if action=='ReLU':
			self.K=np.array([max(0,elem) for elem in self.Kread])
			return

		(imin,imax,idip,izero)=self.IdentifyMaxK()
		if action=='linear':
			interpGap=interp1d([self.times[imax],self.times[izero]], [self.Kread[imax],self.Kread[izero]], kind='linear')
			self.K=np.copy(self.Kread)
			self.K[imax:izero]=interpGap(self.times[imax:izero])
			return

		if action=='exp':
			self.K=self.Kread.copy()
			tau=(self.times[izero]-self.times[imax])/np.log(self.Kread[imax]/self.Kread[izero])
			A=np.exp(self.times[imax]/tau)*self.Kread[imax]
			print('A =',A,'\ttau =',tau, 'Kread[imax] =',self.Kread[imax], 'Kread[izero] =',self.Kread[izero])
			self.K[imax:izero]=A*np.exp(-self.times[imax:izero]/tau)
			return

		if action=='custom':
			self.K=self.Kread.copy()
			imax+=40
			izero+=100
			tau=(self.times[izero]-self.times[imax])/np.log(self.Kread[imax]/self.Kread[izero])
			A=np.exp(self.times[imax]/tau)*self.Kread[imax]
			print('A =',A,'\ttau =',tau, 'Kread[imax] =',self.Kread[imax], 'Kread[izero] =',self.Kread[izero])
			self.K[imax:izero]=A*np.exp(-self.times[imax:izero]/tau)
			return



	def PlotK(self):
		fig, ax = plt.subplots(1,1)
		ax.set_xscale('log')
		ax.set_xlabel('$t$')
		ax.set_ylabel('$\mathcal{K}(t)$')
		ax.plot(self.times, self.K, 'o', label='$\mathcal{K}$')
		ax.plot(self.times, self.Kread, label='$\mathcal{K}_\mathrm{read}$')
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
		ax4.plot(self.times, self.K, color='darkred')
		plt.show()

	def PlotCPP(self):
		fig, ax = plt.subplots(1,1)
		ax.set_ylabel('Correlation Function')
		ax.set_xlabel('$t$')
		ax.set_xscale('log')
		ax.plot(self.times, self.CPP.item()['mean'], color='darkblue' , label='$C^{PP}(t)$')
		ax.plot(self.times, self.CPPdot            , color='darkgreen', label='$\dot{C}^{FP}(t)$')
		ax.plot(self.times, -self.CFP.item()['mean'], color='darkred'  , label='$-C^{FP}(t) = C^{PF}(t)$')
		ax.legend()
		plt.show()


	def CalculateCPPdot(self):
		self.CPPdot=np.ndarray(self.nt, dtype=np.float64)
		f=self.CPP.item()['mean']
		t=self.times
		dfp=np.float64(f[1]-f[0])
		dtp=np.float64(t[1]-t[0])
		self.CPPdot[0]=dfp/dtp
		for i in range(1,self.nt-1):
			dfm=dfp
			dtm=dtp
			dfp=f[i+1]-f[i]
			dtp=t[i+1]-t[i]
			self.CPPdot[i]=0.5*(dfp/dtp+dfm/dtm)
		self.CPPdot[self.nt-1] = dfm/dtm


	def CalculateCFP(self):
		'''
		The noise correlation K is the memory kernel of CPF=-CFP.
		This function calculates CFP using K, in order to compare it to the measured CFP.
		'''
		t=self.times
		cpp=self.CPP.item()['mean']
		interpK=interp1d(t, self.K, kind='cubic')

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
			tempCFP[i]=trapeze(i)*self.invT # The memory kernel is K(t)/kT
		print('')
		self.CFPcheck=tempCFP

	def IntegrateCFP(self):
		self.CPPcheck=np.full(self.nt, self.CPP.item()['mean'][0])
		for i in range(self.nt):
			self.CPPcheck[i]-=integrate.trapz(self.CFPcheck[:i+1], x=self.times[:i+1])
		self.regularizerCPP=self.CPP.item()['mean']/self.CPPcheck

	def PlotCFPcheck(self):
		fig,ax=plt.subplots(1,1)
		ax.set_xscale('log')
		ax.plot(self.times, self.CFPcheck,           color='red'    , label='$C^{FP}$ From memory kernel')
		ax.plot(self.times, self.CFP.item()['mean'], color='darkred', label='$C^{FP}$ Direct measurement')
		ax.plot(self.times, self.K*np.max(self.CFPcheck)/self.K[0], color='grey', label='$\mathcal{K}$')
		ax.plot(self.times, self.CPP.item()['mean']*np.max(self.CFPcheck)/self.CPP.item()['mean'][0], color='lightgrey', label='$C^{PP}$')
		# if hasattr(self,'CFPcheckregularized'):
		# 	ax.plot(self.times, self.CFPcheckregularized,           color='green'    , label='$C^{FP}$ regularized')

		ax.set_title('$T = $%g'%(self.args.temperature))
		ax.set_ylabel('$C(t)$')
		ax.set_xlabel('$t$')
		ax.legend()
		plt.show()

	def PlotCPPcheck(self):
		fig,ax=plt.subplots(1,1)
		ax.set_xscale('log')
		ax.plot(self.times, self.CPPcheck*self.regularizerCPP, 'o',       color='cyan'    , label='$C^{PP}(t) r(t)$ From memory kernel')
		ax.plot(self.times, self.CPPcheck,           color='red'    , label='$C^{PP}(t)$ From memory kernel')
		ax.plot(self.times, self.CPP.item()['mean'], color='darkred', label='$C^{PP}(t)$ Direct measurement')
		ax.set_title('$T = $%g'%(self.args.temperature))
		ax.set_ylabel('$C(t)$')
		ax.set_xlabel('$t$')
		ax.set_ylim((-1,self.CPP.item()['mean'].max()))
		ax.legend()
		plt.savefig('CPPcheckT'+str(self.args.temperature)+'_'+self.args.thermostat+'.png')
		plt.show()

		#Calculate corresponding diffusion coefficient
		istar=np.where(self.times>self.args.tstar)[0][0]
		D=np.trapz(self.CPP.item()['mean'][:istar], x=self.times[:istar])
		Dcheck=np.trapz(self.CPPcheck[:istar], x=self.times[:istar])
		
		nblo=len(self.CPP.item()['blocksJK'])
		Dblocks=np.ndarray(nblo,dtype=np.float64)
		for iblo in range(nblo):
			Dblocks[iblo]=np.trapz(self.CPP.item()['blocksJK'][iblo][:istar], x=self.times[:istar])
		errD=np.sqrt((nblo-1)*(np.square(Dblocks).mean() - D*D  ))

		print('errD = ',errD)

		np.savetxt('CPPcheckT'+str(self.args.temperature)+'_'+self.args.thermostat+'.txt',
			np.column_stack((
				self.times,
				self.CPP.item()['mean'],
				self.CPPcheck, 
				np.full(self.nt, D), 
				np.full(self.nt, errD), 
				np.full(self.nt, Dcheck), 
				np.full(self.nt, self.args.temperature))), 
			header='time CPP CPPcheck D err Dcheck T', fmt="%g %.10g %.10g %.10g %10g %.10g %g")


	def CalculateCFPregularized(self):
		'''
		The noise correlation K is the memory kernel of CPF=-CFP.
		This function calculates CFP using K, in order to compare it to the measured CFP.
		'''
		t=self.times
		cpp=self.CPP.item()['mean']
		interpK=interp1d(t, self.K, kind='cubic')

		def trapeze(i):
			temp=0
			for j in range(1, i+1):
				# print('i:%d  j:%d  ti-tj: %f  ti-tjm1: %f'%(i,j,t[i]-t[j],t[i]-t[j-1]))
				temp += ( interpK(t[i]-t[j])*self.regularizerCPP[j]*cpp[j] + interpK(t[i]-t[j-1])*self.regularizerCPP[j-1]*cpp[j-1] )*(t[j]-t[j-1])
			return 0.5*temp 

		tempCFP=np.ndarray(self.nt ,dtype=np.float64)
		tempCFP[0]=0
		for i in range(0,self.nt):
			print('\rtrapeze iteration', i, end='')
			tempCFP[i]=trapeze(i)*self.invT # The memory kernel is K(t)/kT
		print('')
		self.CFPcheckregularized=tempCFP



if __name__ == "__main__":
	checks=CorrelationConsistency()
	checks.ReadArgs()
	checks.PrintArgs()
	checks.ReadCorrelations()
	# checks.PlotCorrelations()
	checks.CalculateCPPdot()
	# checks.PlotCPP()
	checks.CalculateCFP()
	checks.IntegrateCFP()
	checks.PlotCPPcheck()
	# checks.CalculateCFPregularized()
	# checks.PlotCFPcheck()












