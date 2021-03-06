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
from os import remove,mkdir #Used to remove the backup file at the end of the run
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
from lib import module_measurements as med
from lib import module_timelists as tl
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #


print("hoomd version ",hoomd.__version__)
print("gsd   version ",gsd.__version__)

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
parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='number of atoms in the system')
parser.add_argument('-s','--seed', type=int, required=False, default=-1, help='seed for random numbers. If negative we get it from /dev/urandom')
parser.add_argument('-t','--nNVTsteps', type=int, required=False, default=10000, help='It is the total length of the run (mind that the run often starts at t>0, since we read a backup file). Cannot be negative.')
parser.add_argument('-T','--temperature', type=float, required=False, default=5, help='target Temperature')
parser.add_argument('--tauT', type=float, required=False, default=1.0, help='tau of the thermostat')
parser.add_argument('--dt', type=float, required=False, default=0.002, help='dt for MD integration')
parser.add_argument('--thermostat', required=False, default='NVE', help='basename of the output files')
parser.add_argument('--heavyTrajFreq', type=int, required=False, default=0, help='interval between heavy trajectory backups (default:0, means no backups)')
parser.add_argument('--backupFreq',type=int, required=False, default=0, help='interval between backups (default:0, means no backups)')
parser.add_argument('--trajFreq', type=int, required=False, default=0, help='save trajectory every trajFreq steps (default:0, means no trajectory). If negative, use a logarithmic succession of times, where -trajFreq is the number of configurations in the trajectory (or slightly less, since some times two logarithmic times correspond to the same integer time step).')
parser.add_argument('--iframe', type=int, required=False, default=0, help='Specify from which frame of the gsd file we want the starting configuration (default:0, the first frame)')
parser.add_argument('-l','--label', required=False, default='thermalized', help='basename for the output files')
parser.add_argument('--pot_mode', required=False, default='xplor', choices=['no_shift','shift','xplor'], help='mode for potential')
parser.add_argument('--pot_type', required=False, default='KA', choices=['KA','KAshort','LJmono'], help='type of potential')
parser.add_argument('--addsteps', action='store_true', help='If activated, nNVTsteps are done from the input configuration. If False, we substract ini_step. [Default: False]')
parser.add_argument('--startfromzero', action='store_true', help='If activated, the time step is set to zero once the configuration is read. [Default: False]')
parser.add_argument('--dumpacc', action='store_true', help='If activated, dumps positions, velocities and accelerations in separate .npy files at trajfreq frequency. [Default: False]')
args = parser.parse_args(more_arguments)


if args.seed>0:
	np.random.seed(args.seed)
print("------------------------------")
print("| Input configuration: ",args.filename)
print("| Natoms             = ",args.Natoms)
print("| seed               = ",args.seed)
print("| nNVTsteps          = ",args.nNVTsteps)
print("| T                  = ",args.temperature)
print("| tauT               = ",args.tauT)
print("| dt                 = ",args.dt)
print("| thermostat         = ",args.thermostat)
print("| label              = ",args.label)
print("| pot_mode           = ",args.pot_mode)
print("| heavyTrajFreq      = ",args.heavyTrajFreq)
print("| backupFreq         = ",args.backupFreq)
print("| trajFreq           = ",args.trajFreq)
print("| startfromzero      = ",args.startfromzero)
print("| addsteps           = ",args.addsteps)
print("| dumpacc            = ",args.dumpacc)
print("------------------------------")
if args.nNVTsteps<0:
	raise ValueError('Cannot set a negative number of steps (it is {})'.format(args.nNVTsteps))
if args.dt<0:
	raise ValueError('Cannot set a negative dt (it is {}})'.format(args.dt))
if args.dt>0.05:
	raise ValueError('This dt is too large (it is {}) and the simulations would probably explode'.format(args.dt))
if not args.thermostat=='NVE':
	if args.temperature<0:
		raise ValueError('Cannot set a negative temperature (it is {})'.format(args.temperature))
	if args.tauT<0:
		raise ValueError('Cannot set a negative tauT (it is {})'.format(args.tauT))

################################################################
#
# READ CONFIGURATION
#
################################################################
backupname=args.label+"_backup.gsd"
system = hoomd.init.read_gsd(filename=args.filename, restart=backupname, frame=args.iframe, time_step=(0 if args.startfromzero else None) )
print("| The read configuration has ",len(system.particles)," particles")
Natoms=len(system.particles)
if Natoms!=args.Natoms:
	raise ValueError('The number of particles in the gsd file (Natoms={}) does not match the one given in command line (args.Natoms={})'.format(Natoms,args.Natoms))
iniStep=hoomd.get_step()
print("| iframe: ",args.iframe)
print("| Initial step: ",iniStep)
print("------------------------------")


################################################################
# 
# SET UP POTENTIAL
#
################################################################
print("| Setting potential. type=%s , mode=%s"%(args.pot_type,args.pot_mode))
NeighborsListLJ = md.nlist.cell()
potential=pot.LJ(NeighborsListLJ, type=args.pot_type, mode=args.pot_mode, r_buff=0, Lm=system.box.Lx*0.5-1e-10)

################################################################
# 
# SET UP ANALYZER
#
################################################################
print("| Setting up Analyzer...",end='\t')

#Name of the log
logname=args.label+".txt"

#These are the observables we want to log
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']

print("| In ",logname," we write ",analyzerManyVariables_quantities)

#Every how many integrations we want to log the observables
analyzer_period=int(5./args.dt) #Take measurements once every 5 Lennard Jones times
analyzerManyVariables = hoomd.analyze.log(filename=logname, \
										  quantities=analyzerManyVariables_quantities, period=analyzer_period, \
										  header_prefix = '#seed:'+str(args.seed)+"\n#", \
										  overwrite=False,
										  phase=0)


################################################################
# 
# INTEGRATION
# 
################################################################

def SavePosVelAcc():
	snap=system.take_snapshot()
	np.save(outPos, np.array(snap.particles.position))
	np.save(outVel, np.array(snap.particles.velocity))
	np.save(outAcc, np.array(snap.particles.acceleration))
	return

if args.dumpacc:
	namePos='./pos'+args.label+'.npy'
	nameVel='./vel'+args.label+'.npy'
	nameAcc='./acc'+args.label+'.npy'
	outPos=open(namePos,'ab')
	outVel=open(nameVel,'ab')
	outAcc=open(nameAcc,'ab')


runSteps = max(0,args.nNVTsteps-iniStep) if args.addsteps==False else args.nNVTsteps #If negative, we run no steps
print("| runSteps = ",runSteps)


def SaveGsdNpy(timestep):
	'''
	Saves a gsd file with the current state, and npy files with positions(redundant), velocities(redundant) and Accelerations.
	'''
	SavePosVelAcc()
	hoomd.dump.gsd(filename='trajectory'+args.label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(), phase=-1)
	return

md.integrate.mode_standard(dt=args.dt)
md.update.zero_momentum(phase=-1)
if args.backupFreq>0:
	hoomd.dump.gsd(filename=backupname, overwrite=True, truncate=True, period=args.backupFreq, group=hoomd.group.all(), phase=0)
if args.heavyTrajFreq>0:
	hoomd.dump.gsd(filename='heavyTraj.gsd', overwrite=False, period=args.heavyTrajFreq, group=hoomd.group.all())
if args.trajFreq>0:
	hoomd.dump.gsd(filename='trajectory'+args.label+'.gsd', overwrite=False, period=args.trajFreq, group=hoomd.group.all(),phase=0)
elif args.trajFreq<0:
	nt=-args.trajFreq
	it=0
	listat=tl.ListaLogaritmica(1, runSteps-1, nt, ints=True, addzero=True)
	nt=len(listat) #Since it's a logarithmic list of integers, it might end up having less elements than declare

	#add a trailing element to the list, to avoid IndexError when the list is reached
	listat=np.append(listat,[runSteps])
	if  args.dumpacc:
		hoomd.analyze.callback(callback=SaveGsdNpy, period = lambda n: listat[n] )
	else:
		#If we don't need the accelerations, everything is easy
		hoomd.dump.gsd(filename='trajectory'+args.label+'.gsd', overwrite=True, period = lambda n: listat[n], group=hoomd.group.all() )
	print('| nt=',nt)
	# print("| listat:", listat)



if args.thermostat == 'NVT' :
	print("| ",runSteps,"NVT steps with the Nose-Hoover thermostat at T=",args.temperature)
	integrator = md.integrate.nvt(group=hoomd.group.all(), kT=args.temperature, tau=args.tauT)
	md.update.zero_momentum(period=10./args.dt, phase=0)
	hoomd.run(runSteps, quiet=False)

elif args.thermostat == 'NVE' :
	print("| ",runSteps,"NVE steps with the NVE thermostat")
	integrator = md.integrate.nve(group=hoomd.group.all())
	md.update.zero_momentum(period=10,phase=0)
	hoomd.run(runSteps, quiet=False)

elif args.thermostat == 'MB' :
	if args.trajFreq>=0:
		print("| ",runSteps,"NVT steps with the Andersen thermostat at T=",args.temperature)
		stepsTauT = int(args.tauT/args.dt)
		md.update.zero_momentum(period=10,phase=0)
		integrator = md.integrate.nve(group=hoomd.group.all())
		while(hoomd.get_step()<runSteps):
			for iterations in range(0,int(args.nNVTsteps/stepsTauT)):
				snap = system.take_snapshot()
				vel = np.random.normal(0,np.sqrt(args.temperature), (Natoms,3)) #each component is a gaussian of variance sigma
				vel *= np.sqrt(args.temperature)/np.std(vel,axis=0)
				vel-=vel.mean(axis=0)
				snap.particles.velocity[:] = vel
				system.restore_snapshot(snap)
				hoomd.run(stepsTauT, quiet=False)
	else:
		print("Logarithmic times (trajFreq<0) is not implemented for the MB (Andersen) thermostat")

else:
   print("thermostat=",args.thermostat," is not implemented. EXIT.")
   raise SystemExit

############
# Finalize #
############
#Disable integrator
integrator.disable()

#Save final configuration
finalstatename=args.label+".gsd"
print("finalstatename=",finalstatename)
#Write final state
finalstep=0 if args.thermostat=='MB'else None
hoomd.dump.gsd(filename=finalstatename, overwrite=True, truncate=True, period=None, group=hoomd.group.all(),time_step=finalstep)

#Remove backup because:
# -we don't need to waste memory
# -if the final state is there I need further simulations to start from the final state, not the backup
try:
	remove(backupname)
except OSError: #If the file does not exist
	print("No backup has been generated in this run.")
	pass

#Close I/O streams
if args.dumpacc:
	outPos.close()
	outVel.close()
	outAcc.close()

print("Finished!\n\n")
