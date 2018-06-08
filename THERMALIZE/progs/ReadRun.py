#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This program reads a configuration and runs the MD for a while.
#
################################################################
from __future__ import print_function #for compatibility between python 2 and 3
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys
import os 
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
from lib import module_measurements as med
from lib import module_timelists as tl
from lib import module_signals as sig
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #

################################################################
#
# FUNCTIONS AND CLASSES
# 
################################################################

#HOOMDMAXSTEPS=4000000100
HOOMDMAXSTEPS=4000000000 #Hoomd does not support more than 2^32-1 steps, which we truncate to 4x10^9 for eye friendliness

from pandas import read_csv



#------------------------------------------------#

class Cparams:
	"""
	Class containing all the parameters of the simulation.
	Constructor takes as input the command line arguments.
	Most of the parameters are contained in a file whose name is given through command line.
	"""
	def __init__(self,arguments):
		#Start with command line arguments
		parser = argparse.ArgumentParser(add_help=True)
		parser.add_argument('initconfname', type=str, help='name of the .gsd configuration we want to read')
		parser.add_argument('-p','--paramfilename', type=str, required=False, default=None, help='name of the file with the simulation parameters')
		parser.add_argument('-l','--label', required=False, type=str, default='thermalized', help='basename for the output files')
		parser.add_argument('--iframe', type=int, required=False, default=0, help='From which frame of the gsd file we want the starting configuration (default:0, the first frame)')
		args = parser.parse_args(arguments)
		self.initconfname = args.initconfname
		self.paramfilename = args.paramfilename
		self.label = args.label
		self.iframe = args.iframe

		#Read the file containing the simulation parameters
		try:
			params=np.array(read_csv(self.paramfilename, sep=' ', comment='#',header=None))
		except (ValueError, AttributeError):
			print("No parameters filename given. Running with default bullshit parameters")
		else:
			paramdict = {rows[0]:rows[1] for rows in params}
		self.filename = self.paramfilename
		self.Natoms = 65         if self.paramfilename==None else np.int32(paramdict['Natoms']) #Number of particles
		self.seed = 12345        if self.paramfilename==None else np.int32(paramdict['seed']) #seed for random numbers. If negative we get it from /dev/urandom
		self.nSteps = 1000       if self.paramfilename==None else np.int64(paramdict['nSteps']) #Total length of the run
		self.temperature = 5.0   if self.paramfilename==None else np.float(paramdict['temperature']) 
		self.dt=0.0025           if self.paramfilename==None else np.float(paramdict['dt'])
		self.thermostat='NVT'    if self.paramfilename==None else str(paramdict['thermostat'])
		self.tauT=0.1            if self.paramfilename==None else np.float(paramdict['tauT']) #Tau of the thermostat
		self.backupFreq=0        if self.paramfilename==None else int(np.float64(paramdict['backupFreq'])) #interval between backups (default:0, means no backups)
		self.trajFreq=0          if self.paramfilename==None else int(np.float64(paramdict['trajFreq'])) #save trajectory every trajFreq steps (default:0, means no trajectory). If negative, use a logarithmic succession of times, where -trajFreq is the number of configurations in the trajectory (or slightly less, since some times two logarithmic times correspond to the same integer time step)
		self.addsteps=False      if self.paramfilename==None else (False if paramdict['addsteps']     =='False' else True) #If True, nNVTsteps are done from the input configuration. If False, we substract ini_step. [Default: False]
		self.startfromzero=False if self.paramfilename==None else (False if paramdict['startfromzero']=='False' else True) #If False, initial step is read from the gsd configuration. If True, it is set to zero. [Default: False]
		self.endatzero=False     if self.paramfilename==None else (False if paramdict['endatzero']    =='False' else True) #If True, the time_step of the final configuration is set to zero. [Default: False]
		self.potential='KAshort' if self.paramfilename==None else str(paramdict['potential']) #If True, the time_step of the final configuration is set to zero. [Default: False]
		heavyTrajFreq=0          if self.paramfilename==None else np.int64(np.float64(paramdict['heavyTrajFreq'])) #interval between heavy trajectory backups (default:0, means no backups)
		dumpAcc=True             if self.paramfilename==None else (False if paramdict['dumpAcc']    =='False' else True) #If True, heavy trajectory dumps accelerations. [Default: True]
		maxtimesHT=40            if self.paramfilename==None else np.int32(paramdict['maxtimesHT']) #Number of dumps in the heavy trajectory
		self.heavyTraj=HeavyTraj(freq=heavyTrajFreq, dumpAcc=dumpAcc,maxtimes=maxtimesHT)

		#Derived or hardcoded parameters
		self.logname=self.label+".txt"
		self.analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
		self.analyzer_period = int(10./self.dt) #Take measurements once every 5 Lennard Jones times
		self.backupName=self.label+"_backup.gsd"
		self.heavyTrajDir='./heavyTraj/'
		self.trajName=self.label+'_trajectory.gsd'
		self.timestepName=self.label+'_timestep.in'
		self.finalstatename=self.label+".gsd"
		self.nt= -self.trajFreq if self.trajFreq<0 else np.nan

		self.InitSeed()
		self.Show()
		self.CheckConsistency()

	def CheckConsistency(self):
		if self.Natoms<=0: raise ValueError("Natoms="+str(self.Natoms)+" cannot be negative")
		if not self.nSteps > 0: raise ValueError("nSteps="+str(self.nSteps)+" cannot be negative")
		if not self.temperature > 0: raise ValueError("temperature="+str(self.temperature)+" cannot be negative")
		if not self.dt > 0: raise ValueError("dt="+str(self.dt)+" cannot be negative")
		if not self.dt < 0.01: raise ValueError("dt="+str(self.dt)+" is too big")
		if self.thermostat not in ['NVT','NVE','MB']: raise ValueError("thermostat="+str(self.thermostat)+" should be NVT, NVE or MB")
		if not self.tauT > 0: raise ValueError("tauT="+str(self.tauT)+" cannot be negative")
		return

	def InitSeed(self):
		if self.seed<=0:
			self.seed=np.random.randint(2**32)
		np.random.seed(self.seed)
		return

	def Show(self):
		print("#-----------------------------------------------------#")
		print("# Simulation directory : ",os.getcwd())
		print("# Parameter File Name  : ",self.paramfilename)
		print("# Timestep File Name   : ",self.timestepName)
		print("# Backup File Name     : ",self.backupName)
		print("# Trajectory File Name : ",self.trajName)
		print("# Input configuration  : ",self.initconfname)
		print("# iframe               : ",self.iframe)
		print("# Natoms               = ",self.Natoms)
		print("# potential            = ",self.potential)
		print("# seed                 = ",self.seed)
		print("# nSteps               = ",self.nSteps)
		print("# T                    = ",self.temperature)
		print("# tauT                 = ",self.tauT)
		print("# dt                   = ",self.dt)
		print("# thermostat           = ",self.thermostat)
		print("# label                = ",self.label)
		print("# backupFreq           = ",self.backupFreq)
		print("# trajFreq             = ",self.trajFreq)
		print("# addsteps             = ",self.addsteps)
		print("# startfromzero        = ",self.startfromzero)
		self.heavyTraj.Show()
		return


#------------------------------------------------#

class HeavyTraj:
	def __init__(self, directory='./heavyTraj/', freq=0, dumpAcc=True, maxtimes=40):
		if freq>HOOMDMAXSTEPS:
			raise OverflowError('Maximum admitted freq for Heavy Trajectory is %d. Change the settings in the parameters file (usually params.in)'%HOOMDMAXSTEPS)
		if freq<0: raise ValueError('Found a negative heavyTrajFreq:'+str(freq))
		self.freq=np.int64(freq)        # How often I record a heavy trajectory
		if freq>0:
			self.dir=directory    # Output directory for the heavy trajectories
			try: os.mkdir(self.dir)
			except FileExistsError: pass
			self.len=self.freq-1                    # Length of each heavy trajectory
			self.t0=None                            # Starting time of the ongoing Heavy Trajectory
			self.nextt0=None                        # Starting time of the next Heavy Trajectory
			self.nextt=None                         # Next time at which dump is needed
			self.nextit=None                        # Index of next time at which dump is needed
			self.nameTimes=self.dir+'/times.txt'    # Name of the times
			self.namePos=self.dir+'/pos.npy'        # Name of the positions
			self.nameVel=self.dir+'/vel.npy'        # Name of the velocites
			self.nameAcc=self.dir+'/acc.npy'        # Name of the accelerations
			self.maxtimes=maxtimes                  # How many configurations in a heavy trajectory
			self.ntimes=None                        # How many configurations in a heavy trajectory
			self.outTimes=open(self.nameTimes,'at') # Stream for writing the time list
			self.dumpAcc=dumpAcc                    # Whether to dump accelerations besides the positions
			self.outPos=open(self.namePos,'ab')     # Stream for writing the time list
			if dumpAcc:
				self.outVel=open(self.nameVel,'ab') # Stream for writing the time list
				self.outAcc=open(self.nameAcc,'ab') # Stream for writing accelerations

	def Initt0(self, current_step):
		if self.freq > 0:
			i0    = current_step//self.freq
			self.t0        = np.int64(i0*self.freq)
			inext = i0+1
			self.nextt0    = np.int64(inext*self.freq)
			self.timelist=np.int64(self.t0)+tl.ListaLogaritmica(1, self.len, self.maxtimes, ints=True, addzero=True) if self.freq>0 else None
			print('timelist: ',self.timelist)
			self.ntimes=len(self.timelist)
			#Now find at which element of timelist we are currently
			for i in range(self.ntimes):
				if self.timelist[i] >= current_step: break
			self.nextit=i
			self.nextt=self.timelist[self.nextit]
		return

	def Show(self):
		'''Print all the variables of the class'''
		print("#-----------------------------------------------------#")
		print("# - Heavy Trajectory -")
		for attr in vars(self):
			print("# ",attr,": ",getattr(self,attr))
		print("#-----------------------------------------------------#")
		return

#------------------------------------------------#

class Csim:
	def __init__(self):
		self.InitContext()
		self.params=Cparams(self.more_arguments)
		self.InitConf()
		self.InitSteps()
		self.InitPotential()
		self.InitAnalyzer()
		self.InitIntegrator()
		self.InitDumps()

	def InitContext(self):
		self.context=hoomd.context.initialize()
		hoomd.option.set_notice_level(0)
		self.more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd
		return

	def InitConf(self):
		self.params.iniStep=0 if self.params.startfromzero==True else None
		self.system = hoomd.init.read_gsd(filename=self.params.initconfname, restart=self.params.backupName, frame=self.params.iframe, time_step=self.params.iniStep)
		if not self.params.Natoms==len(self.system.particles): raise SystemExit("Read configuration has Natoms="+str(len(self.system.particles))+"!="+str(self.params.Natoms))
		self.params.iniStep = np.int64(hoomd.get_step())
		return

	def InitPotential(self):
		self.NeighborsList = md.nlist.cell()
		if self.params.Natoms<500: assert(self.params.potential=='KAshort')
		self.pair=pot.LJ(self.NeighborsList, type=self.params.potential)
		return

	def InitAnalyzer(self):
		self.analyzer = hoomd.analyze.log(filename=self.params.logname, \
											  quantities=self.params.analyzer_quantities, period=self.params.analyzer_period, \
											  header_prefix = '#seed:'+str(self.params.seed)+"\n#", \
											  overwrite=False,
											  phase=0)
		return

	def InitIntegrator(self):
		self.mode=md.integrate.mode_standard(dt=self.params.dt)
		if   self.params.thermostat == 'NVT': self.integrator = md.integrate.nvt(group=hoomd.group.all(), kT=self.params.temperature, tau=self.params.tauT)
		elif self.params.thermostat == 'NVE': self.integrator = md.integrate.nve(group=hoomd.group.all())
		elif self.params.thermostat == 'MB' : 
			self.integrator = md.integrate.nve(group=hoomd.group.all())
			if self.params.trajFreq<0: raise SystemExit("Logarithmic times (trajFreq<0) is not implemented for the MB (Andersen) thermostat")
		else: raise SystemExit("thermostat=",self.params.thermostat," is not implemented. EXIT.")
		md.update.zero_momentum(phase=-1)
		return

	def InitSteps(self):
		'''
		Function that deals with the amount of steps that need to be done in this run
		'''
		#See how many steps have been performed previously, that are not stored in the gsd
		try: timedata=read_csv(self.params.timestepName, sep=' ', comment='#')
		except FileNotFoundError: 
			killer = sig.GracefulKiller()
			print(self.params.timestepName," does not exist, so we:\n \
				x) create "+self.params.timestepName+"\n\
				x) fill it with the time step we read in the initial gsd configuration\n\
				x) set to zero the time step on the gsd configuration\n\
				x) Exit.\n\
				xx) Now you can rerun and it should proceed smoothly."
				)
			self.params.numFullCycles=0
			self.params.iCycle=0
			self.params.stepsPerCycle=hoomd.get_step()
			self.params.totOldCycleStep=hoomd.get_step() #steps done in previous cycles
			self.CreateTimeDataFile()
			hoomd.dump.gsd(filename=self.params.finalstatename, overwrite=True, truncate=True, period=None, time_step=0, group=hoomd.group.all())
			self.RemoveBackup()
			raise SystemExit
		else: 
			self.numFullCycles=np.int64(len(timedata))
			if self.numFullCycles==0: raise SystemExit(self.params.timestepName+" looks empty (at most it has the header). Please delete it and then rerun.")
			self.params.iCycle=np.int64(timedata['iCycle'][self.numFullCycles-1]) #take the last line of the file
			self.params.stepsPerCycle=np.int64(timedata['stepsPerCycle'][self.numFullCycles-1])
			self.params.totOldCycleStep=np.int64(timedata['totSteps'][self.numFullCycles-1]) #steps done in previous cycles


		#How many steps need to be done, considering the ones already done, and the addsteps option (which tells you to forget about the past steps)
                if self.params.addsteps==False:
			self.params.runSteps = self.params.nSteps % HOOMDMAXSTEPS
                else:
			totStepsRemaining=self.params.nSteps-(self.params.totOldCycleStep+self.params.iniStep)
			self.params.runSteps = totStepsRemaining% HOOMDMAXSTEPS - self.params.iniStep

		print("#-----------------------------------------------------#")
		print("# initial step      = ",self.params.iniStep)
		print("# older steps       = ",self.params.totOldCycleStep)
		print("# total steps simulated for this sample = ", self.params.totOldCycleStep + self.params.iniStep)
		print("# total steps requested for this sample = ", self.params.nSteps)
		print("# steps to do in this run =",self.params.runSteps)
		print("#-----------------------------------------------------#")

		return

	def CreateTimeDataFile(self):
		with open(self.params.timestepName,"wt") as outf:
			print("iCycle stepsPerCycle totSteps",file=outf)
			print("%g %g %g"%(self.params.iCycle,self.params.stepsPerCycle,self.params.totOldCycleStep),file=outf)
		return

	def AppendTimeDataFile(self):
		'''
		Write the final time in the timestep.in file, so that when it is read again we know
		how long the configuration has been thermalized.
		This crap is necessary because shitty hoomd does not handle runs of more than 2^32~10^9 steps (and I need like 10^13)
		'''
		newiCycle        = self.params.iCycle + 1
		newStepsPerCycle = hoomd.get_step()
		totCycleStep     = np.int64(newStepsPerCycle) + self.params.totOldCycleStep
		with open(self.params.timestepName,"at") as outf:
			print("%g %g %g"%(newiCycle,newStepsPerCycle,totCycleStep),file=outf)
		return

	def InitDumps(self):
		'''
		Initialize settings for writing configurations to disk.
		- backup: in case the run is interrupted
		- trajectory: save evolution of a trajectory (can be logarithmically spaced)
		- heavy trajectory:
			-at each time t0 create a new trajectory of length heavyTrajFreq and save it
		'''

		# BACKUPS
		if self.params.backupFreq>0:
			self.backupDump=hoomd.dump.gsd(filename=self.params.backupName, overwrite=True, truncate=True, period=self.params.backupFreq, group=hoomd.group.all(), phase=0)

		# TRAJECTORY
		if self.params.trajFreq>0:
			self.trajDump=hoomd.dump.gsd(filename=self.params.trajName, overwrite=False, period=self.params.trajFreq, group=hoomd.group.all(),phase=0)
		elif self.params.trajFreq<0:
			self.params.listat=tl.ListaLogaritmica(1, self.params.runSteps, self.params.nt, ints=True, addzero=True)
			self.trajDump=hoomd.dump.gsd(filename=self.params.trajName, overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
			self.params.nt=len(self.params.listat) #Since listat is a logarithmic list of integers, it might end up having less elements than declared

		#HEAVY TRAJECTORY
		self.runCallback=None
		self.runCallbackFreq=0
		if self.params.heavyTraj.freq>0:
			if self.params.trajFreq<0: print('WARNING - HEAVY TRAJECTORY IS NOT IMPLEMENTED FOR trajFreq<0')
			if self.params.thermostat == 'MB' : print('WARNING - HEAVY TRAJECTORY IS NOT IMPLEMENTED FOR MB thermostat')
			totaltime = self.params.totOldCycleStep+self.params.iniStep
			self.params.heavyTraj.Initt0(totaltime) #This is needed in case the previous run didn't complete the heavy trajectory
			# self.heavyTrajPhase = (self.params.heavyTraj.nextt0 - totaltime)%self.params.heavyTraj.freq
			self.runCallback=self.MakeHeavyTrajectory
			self.runCallbackFreq=1
			with open(self.params.heavyTraj.dir+'L.txt','w') as outL:
				outL.write("%.16g\n"%self.system.box.Lx)
				outL.write("%.16g\n"%self.params.dt)			
		return

	def MakeHeavyTrajectory(self,timestep):
		totaltime = self.params.totOldCycleStep+np.int64(timestep)

		if totaltime%self.params.heavyTraj.freq == 0:
			self.params.heavyTraj.Initt0(totaltime)
			# np.savetxt(self.params.heavyTraj.outTimes, np.column_stack(( self.params.heavyTraj.t0*np.ones( self.params.heavyTraj.ntimes), self.params.heavyTraj.timelist)),fmt='%lu')

		if totaltime==self.params.heavyTraj.nextt:
			killer = sig.GracefulKiller()

			#save time
			self.params.heavyTraj.outTimes.write("%d %d\n"%(self.params.heavyTraj.t0,totaltime))
			self.params.heavyTraj.outTimes.flush()

			#save positions, velocities and accelerations
			snap=self.system.take_snapshot()
			np.save(self.params.heavyTraj.outPos, np.array(snap.particles.position))
			if self.params.heavyTraj.dumpAcc==True:
				np.save(self.params.heavyTraj.outVel, np.array(snap.particles.velocity))
				np.save(self.params.heavyTraj.outAcc, np.array(snap.particles.acceleration))

			#Update iterators
			self.params.heavyTraj.nextit=self.params.heavyTraj.nextit+1
			if self.params.heavyTraj.nextit < self.params.heavyTraj.ntimes:
				self.params.heavyTraj.nextt = self.params.heavyTraj.timelist[self.params.heavyTraj.nextit]

			if killer.kill_now: raise SystemExit("Termination signal was caught. Exiting Gracefully...")
			killer.resetSignals() #Now that we are past the dangerous zone, we go back to normal signal handling
			del killer
		return

	def Run(self):
		'''
		Almost all thermostats can be run with the same syntax, except MB, that needs some works, so I had
		to define a different function.
		'''
		md.update.zero_momentum(period=int(1./self.params.dt),phase=0)
		if   self.params.thermostat == 'NVT': self.RunNormal()
		elif self.params.thermostat == 'NVE': self.RunNormal()
		elif self.params.thermostat == 'MB' : self.RunMB()
		else: raise SystemExit("thermostat=",self.params.thermostat," is not implemented. EXIT.")
		return

	def RunNormal(self):
		'''
		Function that runs the dynamics.
		If trajFreq>=0, it just does hoomd.run()
		otherwise, it takes care of saving the trajectory in logarithmic steps.
		'''
		print("# ",self.params.runSteps," steps with the ",self.params.thermostat," thermostat at T=",self.params.temperature)
		sys.stdout.flush()
		if self.params.trajFreq>=0:
#			self.params.runSteps=10000
			print(self.params.runSteps,'steps to do')
			hoomd.run(int(self.params.runSteps), quiet=False, callback=self.runCallback, callback_period=self.runCallbackFreq)
#			hoomd.run(10000, quiet=False, callback=self.runCallback, callback_period=self.runCallbackFreq)
			print(self.params.runSteps,'steps done')
		else:
			for it in range(self.params.nt-1):
				curstep = hoomd.get_step()-self.params.iniStep
				if curstep>=self.params.listat[it+1]:
					continue
				fewSteps=self.params.listat[it+1]-curstep;#self.params.listat[it]
				hoomd.run(fewSteps, quiet=False)
				hoomd.dump.gsd(filename=self.params.trajName, overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
				it+=1
		return

	def RunMB(self):
		'''
		Runs dynamics with Andersen thermostat
		'''
		print("# ",self.params.runSteps," MB steps with the Andersen thermostat at T=",self.params.temperature)
		sys.stdout.flush()
		if self.params.trajFreq>=0:
			stepsTauT = min(int(self.params.tauT/self.params.dt),self.params.runSteps)
			for it in range(0,int(self.params.runSteps/stepsTauT)):
				print('it: ',it,'stepsTauT=',stepsTauT)
				self.integrator.randomize_velocities(kT=self.params.temperature, seed=np.random.randint(2**32))
				print('randomize_velocities finished')
				hoomd.run(stepsTauT, quiet=False)
				print('run finished')
		else: 
			raise SystemExit("Logarithmic times (trajFreq<0) is not implemented for the MB (Andersen) thermostat")

	def RemoveBackup(self):
		"""
		Remove backup because:
		-no need to waste memory
		-if the final state is there I need further simulations to start from the final state, not the backup
		"""
		try: os.remove(self.params.backupName)
		except OSError: print("No backup has been generated in this run.")
		return

	def Finalize(self):
		'''
		Do some final stuff for the simulation, such as saving configuration and removing conflicting files.
		'''
		#First, we make sure that this function is not interrupted
		#When initiated, the class GracefulKiller catches all the sigterm,sigint
		#and sets to True the flag kill_now
		killer = sig.GracefulKiller()

		#If the cycle is finished (we did the maximum number of steps allowed by hoomd)
		#then we need to set the count to zero, and update the configuration
		curTimeStep=hoomd.get_step()
		if (curTimeStep + (self.params.iniStep%HOOMDMAXSTEPS) >= HOOMDMAXSTEPS): 
			self.params.endatzero=True
			self.AppendTimeDataFile()
		else: 
			print("# step:",hoomd.get_step(),"  maxsteps:",HOOMDMAXSTEPS)

		#disable integrator
		if self.integrator.enabled: self.integrator.disable()

		#save final configuration and remove backup
		timestep=0 if self.params.endatzero==True else curTimeStep
		hoomd.dump.gsd(filename=self.params.finalstatename, overwrite=True, truncate=True, period=None, time_step=timestep,group=hoomd.group.all())
		if os.path.isfile(self.params.finalstatename): self.RemoveBackup()
		else: raise SystemExit('Error writing '+self.params.finalstatename)

		#Close open streams
		if self.params.heavyTraj.freq>0:
			self.params.heavyTraj.outTimes.close()
			self.params.heavyTraj.outPos.close()
			if self.params.heavyTraj.dumpAcc: 
				self.params.heavyTraj.outVel.close()
				self.params.heavyTraj.outAcc.close()

		#If a sigterm was sent, now we can terminate the process
		if killer.kill_now: raise SystemExit("Termination signal was caught. Exiting Gracefully...")
		killer.resetSignals() #Now that we are past the dangerous zone, we go back to normal signal handling
		del killer

		return





################################################################
#
# HERE STARTS THE PROGRAM
# 
################################################################
from lib import module_debug as dbg
# sys.settrace(dbg.trace)

if __name__ == "__main__":
	sim=Csim()
	sim.Run()
	sim.Finalize()
	print("# Finished!")
