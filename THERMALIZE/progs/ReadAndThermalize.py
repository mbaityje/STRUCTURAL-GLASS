#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file at temperature T_ini,
# and thermalizes it at temperature T_final.
#
#
# To display help:
# python ReadAndThermalize.py --user="-h"
#
# To launch a simulation:
# python ReadAndThermalize.py --user="filename --dt=dt <and more arguments>"
#
# List of arguments:
#                filename
# -N             Natoms
# -s             seed (<0: /dev/urandom)
# -t             if>0:Max time steps, if<0: number of NVT steps between thermalization checks
# -T             temperature
# --tauT         tau of the thermostat
# -d             MD integration step dt
# --thermostat   thermostat
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
from os import remove #Used to remove the backup file at the end of the run
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
from lib import module_measurements as med
from lib import module_timelists as tl
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #

################################################################
#
# SET UP THE SIMULATION CONTEXT
# 
################################################################
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

################################################################
#
# READ ARGUMENTS
# Read command line arguments
#
#                filename
# -N             Natoms
# -s             seed (<0: /dev/urandom)
# -t             if>0:Max time steps, if<0: number of NVT steps between thermalization checks
# -T             temperature
# --tauT         tau of the thermostat
# -d             MD integration step dt
# --thermostat   thermostat
# --deltaHeavyTraj  Every deltaHeavyTraj steps a configuration is saved. This configuration
#                   can be used to calculate the self-intermediate scattering function or as
#                   a starting point in case at some point crystallization was reached.
#
################################################################
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', #positional argument
					nargs=1,
					help='name of the .gsd configuration we want to read'
)
parser.add_argument('-N','--Natoms', #mandatory
					nargs=1,
					type=int,
					required=True,
					help='number of atoms in the system'
)
parser.add_argument('-l','--label', #optional argument
					nargs=1,
					required=False,
					default=['thermalized'],
					help='basename for the output files'
)
parser.add_argument('-s','--seed', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[-1],
					help='seed for random numbers. If negative we get it from /dev/urandom'
)
parser.add_argument('-t','--nNVTsteps', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[10000],
					help='It is the total length of the run (mind that the run often starts at t>0, since we read a backup file). Cannot be negative.'
)
parser.add_argument('-T','--temperature', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[5],
					help='target Temperature'
)
parser.add_argument('--tauT', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[1.0],
					help='tau of the thermostat'
)
parser.add_argument('--dt', #optional argument
					nargs=1,
					type=float,
					required=False,
					default=[0.002],
					help='dt for MD integration'
)
parser.add_argument('--thermostat', #optional argument
					nargs=1,
					required=False,
					default=['NVE'],
					help='basename of the output files'
)

parser.add_argument('--heavyTrajFreq', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='interval between heavy trajectory backups (default:0, means no backups)'
)
parser.add_argument('--backupFreq', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='interval between backups (default:0, means no backups)'
)
parser.add_argument('--iframe', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='Specify from which frame of the gsd file we want the starting configuration (default:0, the first frame)'
)
parser.add_argument('--trajFreq', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='save trajectory every trajFreq steps (default:0, means no trajectory). If negative, use a logarithmic succession\
					 of times, where -trajFreq is the number of configurations in the trajectory (or slightly less, since some times two \
					 logarithmic times correspond to the same integer time step).'
)
parser.add_argument('--addsteps', #optional argument
					nargs=1,
					type=bool,
					required=False,
					default=[False],
					help='If True, nNVTsteps are done from the input configuration. If False, we substract ini_step. [Default: False]'
)
parser.add_argument('--startfromzero', #optional argument
					nargs=1,
					type=bool,
					required=False,
					default=[False],
					help='If False, initial step is read from the gsd configuration. If True, it is set to zero. [Default: False]'
)


args = parser.parse_args(more_arguments)
filename=args.filename[0]
Natoms=args.Natoms[0]
label=args.label[0]
seed=args.seed[0]
heavyTrajFreq=args.heavyTrajFreq[0]
backupFreq=args.backupFreq[0]
trajFreq=args.trajFreq[0]
iframe=args.iframe[0]
addsteps=args.addsteps[0]
startFromZero=args.startfromzero[0]

nNVTsteps=args.nNVTsteps[0]

TemperatureGoal=args.temperature[0]
tauT=args.tauT[0]
dt=args.dt[0]
thermostat=args.thermostat[0]
del parser

if seed>0: np.random.seed(seed)
print("Input configuration: ",filename)
print("Natoms = ",Natoms)
print("seed = ",seed)
print("nNVTsteps = ",nNVTsteps)
print("T = ",TemperatureGoal)
print("tauT = ",tauT)
print("dt = ",dt)
print("thermostat = ",thermostat)
print("label = ",label)
assert(nNVTsteps>0)
assert(TemperatureGoal>0)
assert(tauT>0)
assert(dt>0 and dt<0.1)

################################################################
#
# READ CONFIGURATION
#
################################################################
backupname=label+"_backup.gsd"
iniStep=0 if startFromZero==True else None
system = hoomd.init.read_gsd(filename=filename, restart=backupname, frame=iframe, time_step=iniStep)
print("The read configuration has ",len(system.particles)," particles")
assert(Natoms==len(system.particles))
iniStep=hoomd.get_step()
print("iframe: ",iframe)
print("Initial step: ",iniStep)


################################################################
# 
# SET UP POTENTIAL
#
################################################################
NeighborsListLJ = md.nlist.cell()
print(" *** Setting Kob-Andersen Potential *** ")
if Natoms<500:
	myLjPair=pot.KApotentialShort(NeighborsListLJ)
else:
	myLjPair=pot.KApotential(NeighborsListLJ)

################################################################
# 
# SET UP ANALYZER
#
################################################################
print("\n\n\nSET UP ANALYZER\n")

#Name of the log
logname=label+".txt"

#These are the observables we want to log
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']

print("In ",logname," we write ",analyzerManyVariables_quantities)

#Every how many integrations we want to log the observables
analyzer_period=int(5./dt) #Take measurements once every 5 Lennard Jones times
analyzerManyVariables = hoomd.analyze.log(filename=logname, \
										  quantities=analyzerManyVariables_quantities, period=analyzer_period, \
										  header_prefix = '#seed:'+str(seed)+"\n#", \
										  overwrite=False,
										  phase=0)


################################################################
# 
# INTEGRATION
# 
################################################################

runSteps = max(0,nNVTsteps-iniStep) if addsteps==False else nNVTsteps #If negative, we run no steps
print("runSteps = ",runSteps)


md.integrate.mode_standard(dt=dt)
md.update.zero_momentum(phase=-1)
if backupFreq>0:
	hoomd.dump.gsd(filename=backupname, overwrite=True, truncate=True, period=backupFreq, group=hoomd.group.all(), phase=0)
if heavyTrajFreq>0:
	hoomd.dump.gsd(filename='heavyTraj.gsd', overwrite=False, period=heavyTrajFreq, group=hoomd.group.all())
if trajFreq>0:
	hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=trajFreq, group=hoomd.group.all(),phase=0)
elif trajFreq<0:
	nt=-trajFreq
	it=0
	listat=tl.ListaLogaritmica(1, runSteps, nt, ints=True, addzero=True)
	hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
	print("listat:", listat)
	nt=len(listat) #Since it's a logarithmic list of integers, it might end up having less elements than declare




if thermostat == 'NVT' :
	print(runSteps," NVT steps with the Nose-Hoover thermostat at T=",TemperatureGoal)
	integrator = md.integrate.nvt(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)
	md.update.zero_momentum(period=10,phase=0)
	if trajFreq>=0:
		hoomd.run(runSteps, quiet=False)
	else:
		for it in range(nt-1):
			fewSteps=listat[it+1]-listat[it]
			hoomd.run(fewSteps, quiet=False)
			hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
			it+=1

elif thermostat == 'NVE' :
	print(runSteps," NVE steps with the NVE thermostat")
	integrator = md.integrate.nve(group=hoomd.group.all())
	md.update.zero_momentum(period=10,phase=0)
	if trajFreq>=0:
		hoomd.run(runSteps, quiet=False)
	else:
		for it in range(nt-1):
			fewSteps=listat[it+1]-listat[it]
			hoomd.run(fewSteps, quiet=False)
			hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
			it+=1


elif thermostat == 'MB' :
	if trajFreq>=0:
		print(runSteps," NVT steps with the Andersen thermostat at T=",TemperatureGoal)
		stepsTauT = int(tauT/dt)
		md.update.zero_momentum(period=10,phase=0)
		integrator = md.integrate.nve(group=hoomd.group.all())
		while(hoomd.get_step()<runSteps):
			for iterations in range(0,int(nNVTsteps/stepsTauT)):
				snap = system.take_snapshot()
				vel = np.random.normal(0,np.sqrt(TemperatureGoal), (Natoms,3)) #each component is a gaussian of variance sigma
				vel *= np.sqrt(TemperatureGoal)/np.std(vel,axis=0)
				vel-=vel.mean(axis=0)
				snap.particles.velocity[:] = vel
				system.restore_snapshot(snap)
				hoomd.run(stepsTauT, quiet=False)
	else:
		print("Logarithmic times (trajFreq<0) is not implemented for the MB (Andersen) thermostat")

else:
   print("thermostat=",thermostat," is not implemented. EXIT.")
   raise SystemExit

############
# Finalize #
############
integrator.disable()
finalstatename=label+".gsd"
print("finalstatename=",finalstatename)
#Write final state
hoomd.dump.gsd(filename=finalstatename, overwrite=True, truncate=True, period=None, group=hoomd.group.all())
#Remove backup because:
# -we don't need to waste memory
# -if the final state is there I need further simulations to start from the final state, not the backup
try:
	remove(backupname)
except OSError: #If the file does not exist
	print("No backup has been generated in this run.")
	pass
#Finally, if the simulation was run with MB, i want to set the number of time steps to zero,
#because with MB I only do the initial thermalization to a very high temperature.
#Since hoomd does not allow to reset the step number, I need to do a little walkaround.
if thermostat == 'MB':
	hoomd.context.initialize() #I need to create another simulation context, otherwise this shit bitches when I read_gsd
	hoomd.init.read_gsd(filename=finalstatename, restart=None, time_step=0)
hoomd.dump.gsd(filename=finalstatename, overwrite=True, period=None, group=hoomd.group.all())



print("Finished!\n\n")
