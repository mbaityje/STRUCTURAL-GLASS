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
		self.ReadArgs()
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
		parser.add_argument('--smoothK', required=False, default=None, choices=['None','ReLU','linear','exp','zero','custom'], help='Apply smoothing to the K(t)')
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
		self.MSD = np.ndarray.tolist(np.load('msd_NVT.npy'))
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

		# Then find when it touches zero for the first time
		for i in range(1,self.nt):
			if self.Kread[i]<0:
				izero1=i
				break
		#Then find when it cuts zero from the dip
		for i in range(idip+1,self.nt):
			if self.Kread[i]>0:
				izero=i
				break

		print('times[imin] =', self.times[imin],'\ttimes[imax] =', self.times[imax], '\ttimes[idip] =', self.times[idip],'\ttimes[izero] =', self.times[izero],'\ttimes[izero1] =', self.times[izero1])
		return imin,imax,idip,izero, izero1

	def SmoothK(self, action=None):

		print('action = ',action)
		if action==None:
			self.K=self.Kread.copy()
			return

		if action=='ReLU':
			print('Setting to zero all negative elements of K(t)')
			self.K=np.array([max(0,elem) for elem in self.Kread])
			return

		(imin,imax,idip,izero,izero1)=self.IdentifyMaxK() #izero1: when K(t) crosses zero the first time (from positive to negative). izero: when K(t) becomes positive again

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

		if action=='zero':
			print('Setting K(t) to zero as soon as it touches zero')
			self.K=self.Kread.copy()
			self.K[izero1:]=0
			print('izero1 = ',izero1)

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
		fig, (ax1,ax2) = plt.subplots(2,1)
		ax1.set_ylabel('$\Delta^2(t)$')
		ax1.get_xaxis().set_visible(False)
		ax1.set_xscale('log')
		ax1.set_yscale('log')
		ax1.plot(self.times, self.MSD['mean'])
		ax2.set_ylabel('$\mathcal{K}(t)$')
		ax2.set_xlabel('$t$')
		ax2.get_xaxis().set_visible(True)
		ax2.set_xscale('log')
		ax2.plot(self.times, self.K, color='darkred')
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


	def CalculateMSDdot(self):
		self.MSDdot=np.ndarray(self.nt, dtype=np.float64)
		f=self.MSD['mean']
		t=self.times
		dfp=np.float64(f[1]-f[0])
		dtp=np.float64(t[1]-t[0])
		self.MSDdot[0]=dfp/dtp
		for i in range(1,self.nt-1):
			dfm=dfp
			dtm=dtp
			dfp=f[i+1]-f[i]
			dtp=t[i+1]-t[i]
			self.MSDdot[i]=0.5*(dfp/dtp+dfm/dtm)
		self.MSDdot[self.nt-1] = dfm/dtm

	def CalculateMSDdotdot(self):
		self.MSDdotdot=np.ndarray(self.nt, dtype=np.float64)
		f=self.MSDdot
		t=self.times
		dfp=np.float64(f[1]-f[0])
		dtp=np.float64(t[1]-t[0])
		self.MSDdotdot[0]=dfp/dtp
		for i in range(1,self.nt-1):
			dfm=dfp
			dtm=dtp
			dfp=f[i+1]-f[i]
			dtp=t[i+1]-t[i]
			self.MSDdotdot[i]=0.5*(dfp/dtp+dfm/dtm)
		self.MSDdotdot[self.nt-1] = dfm/dtm


	def CalculateMSDdotdotCheck(self):
		'''
		Calculates the second derivative of the MSD through K(t), by using the relation

		MSDdotdotCheck = 6kT/m - int_0^t dt' K(t-t') MSDdot(t')

		The obtained function can be compared with that coming by derivating the measured MSD
		'''
		t=self.times
		UU=self.MSDdot
		KK=self.K*self.invT # The memory kernel is K(t)/kT

		interpK=interp1d(t, KK, kind='cubic')
		interpU=interp1d(t, UU, kind='cubic')

		def trapeze(i):
			temp=0
			for j in range(1, i+1):
				# temp += ( interpK(t[i]-t[j])*UU[j] + interpK(t[i]-t[j-1])*UU[j-1] )*(t[j]-t[j-1])
				
				temp += ( interpU(t[i]-t[j])*KK[j] + interpU(t[i]-t[j-1])*KK[j-1] )*(t[j]-t[j-1])
			return 0.5*temp 


		sixTonM = 6*self.args.temperature #Mass=1, k_B=1
		tempMSDdotdot=np.ndarray(self.nt ,dtype=np.float64)
		tempMSDdotdot[0]=sixTonM
		for i in range(0,self.nt):
			print('\rtrapeze iteration', i, end='')
			tempMSDdotdot[i]=sixTonM - trapeze(i)
		print('')

		self.MSDdotdotCheck=tempMSDdotdot

	def TruncateMSDdotdotCheck(self):

		for i in range(self.nt-1,0,-1):
			if self.MSDdotdotCheck[i]<0:
				break

		itrunc=i#-150
		print('Truncating at itrunc=',itrunc)

		self.MSDdotdotCheck[itrunc:]=0



	def PlotMSDdotdot(self):
		fig,(ax1,ax2)=plt.subplots(2,1)
		ax1.set_xscale('log')
		ax1.set_title('$T = $%g'%(self.args.temperature))
		ax1.plot(self.times, self.MSDdotCheck,           color='red'    , label='From memory kernel')
		ax1.set_ylabel(r'$\frac{d \Delta^2(t)}{dt}$')
		ax1.set_xlabel('$t$')
		ax1.legend()

		ax2.set_xscale('log')
		ax2.plot(self.times, self.MSDdotdotCheck, color='red'    , label='From memory kernel')
		ax2.plot(self.times, self.MSDdotdot, 	  color='darkred', label='Direct measurement')
		ax2.set_ylabel(r'$\frac{d^2 \Delta^2(t)}{dt^2}$')
		ax2.set_xlabel('$t$')
		ax2.legend()
		plt.show()


	def IntegrateMSDdotdot(self):
		self.MSDdotCheck=np.full(self.nt, self.MSDdot[0])
		for i in range(self.nt):
			self.MSDdotCheck[i]=integrate.trapz(self.MSDdotdotCheck[:i+1], x=self.times[:i+1])

		self.MSDCheck=np.full(self.nt, self.MSD['mean'][0])
		for i in range(self.nt):
			self.MSDCheck[i]=integrate.trapz(self.MSDdotCheck[:i+1], x=self.times[:i+1])


	def PlotMSD(self):
		fig,ax=plt.subplots(1,1)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.plot(self.times, self.MSDCheck,           color='red'    , label='From memory kernel')
		ax.errorbar(self.times, self.MSD['mean'], yerr=self.MSD['err'], color='darkred', label='Direct measurement')
		ax.plot(self.times,0.1*self.times, color='grey', label='$t$',linestyle='--')
		ax.plot(self.times,np.square(self.times), color='black', label='$t^2$', linestyle=':')
		ax.set_title('$T = $%g'%(self.args.temperature))
		ax.set_ylabel(r'$\Delta^2(t)$')
		ax.set_xlabel('$t$')
		ax.legend()
		plt.show()


	def FakeMSDdot(self,t, tau=1):

		fourtausq=4*tau*tau
		factor=6*self.args.temperature
		B=np.sqrt(1-1./fourtausq)
		phi = np.arctan2(B,0.5/tau) - np.arctan2(B,-0.5/tau)

		suma=1./tau + np.exp(-0.5*t/tau)*np.sin(B*t+phi)/B

		return factor*suma

	def FakeConsistencyCheck(self):

		#Calculate derivative of MSD from a fake exponential kernel with autocorrelation time tau
		#
		# fakeK(t) = exp(-t/tau)

		tau=2 # Value for T=0.7
		self.fakeMSDdot=[self.FakeMSDdot(self.times[i],tau=tau) for i in range(self.nt)]

		# Find second derivative of fake MSD
		self.fakeMSDdotdot=np.ndarray(self.nt, dtype=np.float64)
		f=self.fakeMSDdot
		t=self.times
		dfp=np.float64(f[1]-f[0])
		dtp=np.float64(t[1]-t[0])
		self.fakeMSDdotdot[0]=dfp/dtp
		for i in range(1,self.nt-1):
			dfm=dfp
			dtm=dtp
			dfp=f[i+1]-f[i]
			dtp=t[i+1]-t[i]
			self.fakeMSDdotdot[i]=0.5*(dfp/dtp+dfm/dtm)
		self.fakeMSDdotdot[self.nt-1] = dfm/dtm


		fig,(ax1,ax2)=plt.subplots(2,1)
		ax1.set_xscale('log')
		ax1.set_xlabel(r'$t$')
		ax1.set_ylabel(r'$\frac{d \Delta^2(t)}{dt}$')
		ax1.plot(self.times, self.fakeMSDdot)
		ax2.set_xscale('log')
		ax2.set_xlabel(r'$t$')
		ax2.set_ylabel(r'$\frac{d^2 \Delta^2(t)}{dt^2}$')
		ax2.plot(self.times, self.fakeMSDdotdot)
		plt.show()


		# Integrate derivative of MSD 
		self.fakeMSD=np.full(self.nt, 0.0)
		for i in range(self.nt):
			self.fakeMSD[i]+=integrate.trapz(self.fakeMSDdot[:i+1], x=self.times[:i+1])




		# Now calculate the second derivative of the MSD through the usual check that involves the
		# memory function. This time with the memory function that we imposed.
		useFakeMSD=True
		t=self.times
		UU=self.fakeMSDdot if useFakeMSD else self.MSDdot
		KK=np.exp(-t/tau) # The memory kernel is K(t)/kT
		interpK=interp1d(t, KK, kind='cubic')
		interpU=interp1d(t, UU, kind='cubic')

		def trapeze(i):
			temp=0
			for j in range(1, i+1):
				# temp += ( interpK(t[i]-t[j])*UU[j] + interpK(t[i]-t[j-1])*UU[j-1] )*(t[j]-t[j-1])
				
				temp += ( interpU(t[i]-t[j])*KK[j] + interpU(t[i]-t[j-1])*KK[j-1] )*(t[j]-t[j-1])
			return 0.5*temp 


		sixTonM = 6*self.args.temperature #Mass=1, k_B=1
		tempMSDdotdot=np.ndarray(self.nt ,dtype=np.float64)
		tempMSDdotdot[0]=sixTonM
		for i in range(0,self.nt):
			print('\rtrapeze iteration', i, end='')
			tempMSDdotdot[i]=sixTonM - trapeze(i)
		print('')
		self.fakeMSDdotdotCheck=tempMSDdotdot


		#Now integrate twice the check on the fake MSD, to obtain the fake MSD
		self.fakeMSDdotCheck=np.full(self.nt, 0.0)
		for i in range(self.nt):
			self.fakeMSDdotCheck[i]=integrate.trapz(self.fakeMSDdotdotCheck[:i+1], x=self.times[:i+1])

		self.fakeMSDCheck=np.full(self.nt, 0.0)
		for i in range(self.nt):
			self.fakeMSDCheck[i]=integrate.trapz(self.fakeMSDdotCheck[:i+1], x=self.times[:i+1])



		fig,ax=plt.subplots(1,1)
		ax.set_xlabel(r'$t$')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylabel(r'$\Delta^2(t)$')
		ax.plot(self.times, self.fakeMSDCheck,'o', label='Fake MSD obtained through fake K(t)')
		ax.plot(self.times, self.fakeMSD, label='Fake MSD')
		ax.errorbar(self.times, self.MSD['mean'], yerr=self.MSD['err'], label='MSD')
		ax.plot(self.times,self.times, color='grey', label='$t$',linestyle='--')
		ax.legend()
		plt.show()




 

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






if __name__ == "__main__":
	checks=CorrelationConsistency()
	checks.PrintArgs()
	checks.ReadCorrelations()
	checks.PlotCorrelations()
	checks.CalculateMSDdot()
	checks.CalculateMSDdotdot()
	checks.CalculateMSDdotdotCheck()


	# checks.TruncateMSDdotdotCheck() # This function truncates the MSDdotdot to zero when it touches zero

	checks.IntegrateMSDdotdot()
	checks.PlotMSDdotdot()
	checks.PlotMSD()

	checks.FakeConsistencyCheck()
	sys.exit()
	# checks.PlotCPP()
	checks.CalculateCFP()
	checks.IntegrateCFP()
	# checks.PlotCPPcheck()
	# checks.CalculateCFPregularized()
	# checks.PlotCFPcheck()












