#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file, and
# let's it evolve with  NVT (first) and NVE (after) dynamics.
#
# To display help:
# python T1-ReadAndEvolve.py --user="-h"
#
# To launch a simulation:
# python T1-ReadAndEvolve.py --user="filename -e nNVEsteps -t nNVTsteps --thermostat thermostat"
#
# For example:
# python T1-ReadAndEvolve.py --user="sample-states/KA_rho=12e-1_N=1e3_T=45e-2_tag=1_type=restart.gsd -e 10 -t 1000 --thermostat NVT"
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from math import sqrt #Just wanna square root
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #

################################################################
#
# FUNCTIONS THAT WOULD GO IN A SEPARATE MODULE
# 
################################################################

def KApotentialShort(NeighborsList):
	eps_AA=1
	eps_AB=1.5
	eps_BB=0.5
	sig_AA=1
	sig_AB=0.8
	sig_BB=0.88
	r_on_cutoff=1.2
	r_cutoff=1.4
	NeighborsList.set_params(r_buff=0.0)
	myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsList)
	myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
	myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
	myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
	myLjPair.set_params(mode="xplor")
	return myLjPair

def KApotential(NeighborsListLJ): #Defines a Kob-Anderson potential
	r_cutoff=2.5
	eps_AA=1 
	eps_AB=1.5
	eps_BB=0.5
	sig_AA=1
	sig_AB=0.8
	sig_BB=0.88
	r_on_cutoff=1.2
	myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
	myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA, r_on=r_on_cutoff*sig_AA)
	myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB, r_on=r_on_cutoff*sig_AB)
	myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB, r_on=r_on_cutoff*sig_BB)
	myLjPair.set_params(mode="xplor")
	return myLjPair





################################################################
#
# SET UP THE SIMULATION CONTEXT
# 
################################################################
print("\n\n\nSET UP THE SIMULATION CONTEXT\n")
# We now describe hoomd's simulation context. In this case we set a single
# context. Look at hoomd documentation for using multiple contexts.
hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
#hoomd.option.set_msg_file(None) #Redirect output to file. None: stdout/stderr 
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd



################################################################
#
# READ ARGUMENTS
# Read command line arguments
# Since this is a tutorial, we use both positional and optional (with flag) arguments:
#
# <filename> is the first argument
# -e             number of NVE steps
# -v             number of NVT steps
# --thermostat   thermostat (NVT or MB)
#
# Note: usually argparse reads from sys.argv. Since hoomd is a piggy
# and wants to have the priority on command line, we had to first read
# the hoomd parameters, and then we use the parser with the remaining
# options, that we saved in more_arguments (and that are input through the
# --user="..." flag).
#
################################################################
print("\n\n\nREAD ARGUMENTS\n")

#The we create a parser for the User files 
parser = argparse.ArgumentParser(prog='python ReadAndEvolve.py [--hoomdi-flags] --user=" HERE YOU PUT THE FOLLOWING ARGUMENTS"',
								 description='The program reads a .gsd configuration and runs it in NVE and NVT',
								 add_help=True)
parser.add_argument('filename', #positional argument
					nargs=1,
					help='name of the .gsd configuration we want to read'
)
parser.add_argument('-e','--nNVEsteps', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='number of NVE steps (default 0)'
)
parser.add_argument('-t','--nNVTsteps', #optional argument
					nargs=1,
					type=int,
					required=False,
					default=[0],
					help='number of NVT steps (default 0)'
)
parser.add_argument('--thermostat', #optional argument
					nargs=1,
					required=False,
					default=['NVT'],
					choices=['NVT', 'MB'],
					help='can be NVT (Nose-Hoover thermostat, default) or MB (Anderson thermostat)'
)
args = parser.parse_args(more_arguments)
filename=args.filename[0]
nNVEsteps=args.nNVEsteps[0]
nNVTsteps=args.nNVTsteps[0]
thermostat=args.thermostat[0]
del parser

print(filename)
print(nNVEsteps)
print(nNVTsteps)
print(thermostat)


################################################################
#
# READ CONFIGURATION
# Read .gsd configuration
#
################################################################
print("\n\n\nREAD CONFIGURATION\n")
# This is the first hoomd method to be invoked (mode=cpu, gpu)
system = hoomd.init.read_gsd(filename=filename)

#We now play a little with the properties of the read configuration
Natoms=len(system.particles)
print("This is the",system.box)
print(system.particles.types)
print("Number of particle types:",len(system.particles.types))
assert(len(system.particles.types)==2),"This tutorial reads binary mixtures"
groupA = hoomd.group.type(name='a-particles', type='A')
groupB = hoomd.group.type(name='b-particles', type='B')
print(len(groupA)," particles of type A")
print(len(groupB)," particles of type B")


id=Natoms-1
print("\nProperties of particle id= ",id)
print("Tag: ",system.particles[id].tag," - Note that id and tag may not coincide (for example if particle 0 is deleted)")
print("Typeid: ",system.particles[id].typeid)
print("Type: ",system.particles[id].type)
print("Position: ",system.particles[id].position)
print("Velocity: ",system.particles[id].velocity)
print("Mass: ",system.particles[id].mass)
print("Diameter: ",system.particles[id].diameter)

print("\nThe properties of the particle can also be printed altogether. Particle 0:")
print(system.particles[0])




################################################################
# 
# SET UP POTENTIAL
#
################################################################
print("\n\n\nSET UP POTENTIAL\n")
#Set up neighbor list
NeighborsListLJ = md.nlist.cell()

#Set Kob-Andersen potential
myLjPair=KApotential(NeighborsListLJ) if Natoms >500 else KApotentialShort(NeighborsListLJ)
#These are two alternative ways of seeing the values we inserted
#print(myLjPair.pair_coeff.get_metadata())
print(myLjPair.pair_coeff.values)

################################################################
# 
# SET UP ANALYZER
#
################################################################
print("\n\n\nSET UP ANALYZER\n")

#Name of the log
logname="test-output/observables-test_"+thermostat+str(nNVTsteps)+"_NVE"+str(nNVEsteps)+".txt"

#These are the observables we want to log
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']

print("In ",logname," we write ",analyzerManyVariables_quantities)

#Every how many integrations we want to log the observables
analyzer_period=100
analyzerManyVariables = hoomd.analyze.log(filename=logname, \
										  quantities=analyzerManyVariables_quantities, period=analyzer_period, \
										  header_prefix = 'This is a test output of the toy program'+sys.argv[0]+'\n', \
										  overwrite=True,
										  phase=0)

 
################################################################
# 
# INTEGRATION
# 
################################################################
print("\n\n\nINTEGRATING DYNAMICS\n")
TemperatureGoal=1
tauT=1
dt=0.001
print("Simulation parameters that we hardcoded in this tutorial:")
print("TemperatureGoal = ",TemperatureGoal)
print("tauT = ",tauT)
print("dt = ",dt)

#There are two integration modes:
# -standard: allows for normal integration methods, such as nve, nvt, MB, brownian, etc...
# -fire: for energy minimization
md.integrate.mode_standard(dt=dt) ##allows Nose-Hoover in NVT, for instance




######## NVT ########
print(nNVTsteps," NVT steps with the ",thermostat," thermostat")
if nNVTsteps > 0 :

	if thermostat == 'NVT' :
		integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)
		hoomd.run(nNVTsteps, quiet=False) 

	elif thermostat == 'MB' :
		stepsTauT = int(tauT/dt)
		integrator_nvt = md.integrate.nve(group=hoomd.group.all())
		for iterations in range(0,int(nNVTsteps/stepsTauT)):
			snap = system.take_snapshot()
			# each component is a gaussian of variance sigma, that's it.
			vel = np.random.normal(0,sqrt(TemperatureGoal), (Natoms,3))
			print('\nVelocities were rescaled by ',TemperatureGoal**0.5/np.std(vel),'\n')
			vel *= sqrt(TemperatureGoal)/np.std(vel) #Riscalo le velocita` per essere sicuro che non ci sono effetti di taglia finita sulla varianza della distribuzione
			snap.particles.velocity[:] = vel
			system.restore_snapshot(snap)
			hoomd.run(stepsTauT, quiet=False)
	else:
		print("In this tutorial the only implemented thermostats are NVT (Nose-Hoover) and MB (Anderson)\n")
		sys.exit()

	integrator_nvt.disable()
		
######## NVE ########
print(nNVEsteps," NVE steps via Velocity-Verlet")
if nNVEsteps > 0 :
	curStep = hoomd.get_step()
	integrator_nve = md.integrate.nve(group=hoomd.group.all())
	hoomd.run(nNVEsteps, profile=False, quiet=False)
	integrator_nve.disable()

analyzerManyVariables.disable()




