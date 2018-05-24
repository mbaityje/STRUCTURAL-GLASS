#!/usr/bin/env python3
#####################################################################
#                                                                   #
# This example shows that nested analyze callbacks cannot be used,  #
# and provides some alternatives.                                   #
#                                                                   #
# Launch as                                                         #
# python T23-Callbacks.py sample-states/N65.gsd                     #
#                                                                   #
#####################################################################

import sys
import hoomd
from hoomd import md
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
if len(sys.argv)==2:
	filename = sys.argv[1]
else:
	print("Launch as:")
	print(sys.argv[0]," configuration.gsd")
hoomd.init.read_gsd(filename=filename)
neighborsList = md.nlist.cell(r_buff=0.0)
pair = md.pair.lj(r_cut=1.4, nlist=neighborsList)
pair.pair_coeff.set('A', 'A', epsilon=1, sigma=1, r_cut=1.4, r_on=1.2)
pair.pair_coeff.set('A', 'B', epsilon=1.5, sigma=0.8, r_cut=1.4*0.8, r_on=1.2*0.8)
pair.pair_coeff.set('B', 'B', epsilon=0.5, sigma=0.88, r_cut=1.4*0.88, r_on=1.2*0.88)
pair.set_params(mode="xplor")
analyzer = hoomd.analyze.log(filename='test.log', quantities=['temperature', 'pressure'], period=100)
mode = md.integrate.mode_standard(dt=0.0025)
integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=2.0, tau=0.1)

class callbacks:
	def __init__(self):
		self.heavyTrajFreq=20
		self.heavyTrajLen=5
		self.backupFreq=10
		self.trajFreq=10
		self.bctime0=0

	def nestedcallbacks(self):
		'''
		This method is to show that nested callbacks give segmentation fault
		'''
		self.heavyTrajCallback = hoomd.analyze.callback(callback = self.localcallback, period = self.heavyTrajFreq)

	def initDumps(self):
		self.backupDump = hoomd.dump.gsd(filename='backup.gsd', overwrite=True, truncate=True, period=self.backupFreq, group=hoomd.group.all())
		self.trajDump = hoomd.dump.gsd(filename='traj.gsd', overwrite=False, period=self.trajFreq, group=hoomd.group.all())

	def localcallback(self,timestep):
		if 'lc' in vars(self):
			print('Repeating callback')
			self.lc.disable()
			#del self.lc
		else:
			print('Running callback for the first time')
		callback_tag='first step: '+str(timestep)
		self.lc=hoomd.analyze.callback(callback = lambda step: print('step ',step,callback_tag), period = 1)

	def branchedcallback(self, timestep):
		if timestep%self.heavyTrajFreq == 0:
			print('New callback')
			if 'bc' in vars(self) and self.bc.enabled: #In case heavyTrajFreq<=heavyTrajLen
				self.bc.disable()
			callback_tag='first step: '+str(timestep)
			self.bc = hoomd.analyze.callback(callback = lambda step: print('step ',step,callback_tag), period = 1)
			self.bctime0 = timestep
		elif timestep==self.heavyTrajLen+self.bctime0:
			self.bc.disable()

	def fcb(self):
		self.bc2 = hoomd.analyze.callback(callback = lambda step: print('step ',step,'first step: ',step // self.heavyTrajLen*self.heavyTrajFreq), period = lambda n: (n % self.heavyTrajLen) + n // self.heavyTrajLen*self.heavyTrajFreq)

nsteps=50
cb=callbacks()
cb.initDumps()

#Set to true to verify that nested callbacks are not supported
youWantSegmentationFault=False
if youWantSegmentationFault:
	cb.nestedcallbacks()
	hoomd.run(nsteps)


# Every cb.heavyTrajFreq time steps save a trajectory of cb.heavyTrajLen steps
# Examples I and II impose cb.heavyTrajFreq = cb.heavyTrajLen
# Examples III and IV have impose cb.heavyTrajFreq > cb.heavyTrajLen
#I) cb.heavyTrajFreq=cb.heavyTrajLen
print('#######')
print('## I ##')
print('#######')
for ht in range(nsteps//cb.heavyTrajFreq):
	cb.localcallback(hoomd.get_step())
	hoomd.run(cb.heavyTrajFreq)

#II) cb.heavyTrajFreq=cb.heavyTrajLen
print('########')
print('## II ##')
print('########')
hoomd.run(nsteps,callback_period=cb.heavyTrajFreq, callback=cb.localcallback)
cb.lc.disable()

#III) cb.heavyTrajFreq>cb.heavyTrajLen
print('#########')
print('## III ##')
print('#########')
for ht in range(nsteps//cb.heavyTrajFreq):
	print('ht:',ht)
	firststep=hoomd.get_step()
	print('New callback')
	lc=hoomd.analyze.callback(callback = lambda step: print('step ',step,'first step: ',firststep), period = 1)
	hoomd.run(cb.heavyTrajLen)
	lc.disable()
	hoomd.run(cb.heavyTrajFreq-cb.heavyTrajLen)

#IV) cb.heavyTrajFreq>cb.heavyTrajLen
print('########')
print('## IV ##')
print('########')
hoomd.run(nsteps,callback_period=1, callback=cb.branchedcallback)

#V) cb.heavyTrajFreq>cb.heavyTrajLen
print('#######')
print('## V ##')
print('#######')
cb.fcb()
hoomd.run(nsteps)



integrator_nvt.disable()


