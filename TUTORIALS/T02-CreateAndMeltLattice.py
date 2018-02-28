#!/usr/bin/python
################################################################
#
#
# DESCRIPTION
# This example creates a lattice configuration and melts it.
# The trajectory is saved so we can watch a cool video.
# I could also measure the crystalline order parameter, to follow the melting,
# but this has not been implemented.
#
# To display help:
# python T2-CreateAndMeltLattice.py --user="-h"
#
# To launch a simulation:
# python T2-CreateAndMeltLattice.py --user="-t nNVTsteps --thermostat thermostat"
#
# For example:
# python T2-CreateAndMeltLattice.py --user="-t 1000 --thermostat NVT"
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
import kobandersen #Here I have put the Kob-Andersen parameters

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
#
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
                                 description='The program reads a .gsd configuration and runs it in NVT',
                                 add_help=True)
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
nNVTsteps=args.nNVTsteps[0]
thermostat=args.thermostat[0]
del parser #I'm new to python, perhaps this is not proper coding

print(nNVTsteps)
print(thermostat)


################################################################
#
# CREATE CONFIGURATION on an LxLxL lattice
#
################################################################
L=6  #Times the unit cell is repeated along each direction
a=1.61 #I choose this value so the rescaling is small (crashes with large rescale)
rho=1.2 #This is a typical working density

print("\n\n\nCREATE CONFIGURATION\n")

#This is the unit cell we define
#We need to define it with 5 elements, because the particle types are assigned at this stage
#For a monodisperse mixture we could use the hoomd built-in unit cells
uc = hoomd.lattice.unitcell(N = 5, #Number of particles in the unit cell
                            a1 = [  a,   0,   0], #Basis vectors of the unit cell
                            a2 = [  0,   a,   0],
                            a3 = [  0,   0,   a],
                            dimensions = 3,
                            position = [[0, 0, 0], [0.5, 0.0, 0.0], [0.0,0.5,0.0], [0.0,0.0,0.5], [0.5,0.5,0.5]],
                            type_name = ['A','A','A','A','B'], #An 80%-20% mixture
                            mass = [1.0, 1.0, 1.0, 1.0, 1.0],
                            charge = [0.0, 0.0, 0.0, 0.0, 0.0],
                            diameter = [1./2, 1./2, 1./2, 1./2, 0.88/2] #I believe that the diameter only plays a role in visualization
);
# lattice made of repeat of unit cell:
system = hoomd.init.create_lattice(unitcell=uc, n=L)
Natoms=len(system.particles)
print("Created a lattice with ",Natoms," particles")
print("This is the",system.box)
print("Volume before rescaling: ",system.box.get_volume())
Box=hoomd.data.boxdim(volume=Natoms/rho)
system.box=Box
print("Volume after rescaling: ",system.box.get_volume())
print("This is the",system.box)

#We dump the initial configuration
hoomd.dump.gsd(filename="./test-output/lattice.gsd", overwrite=True, period=None, group=hoomd.group.all(), time_step=0)
hoomd.dump.gsd(filename="./test-output/trajectory.gsd", overwrite=True, period=50, group=hoomd.group.all(), phase=0)

################################################################
# 
# CREATE GROUPS OF PARTICLES
#
################################################################
#80% of the particles should belong to group A, 20% of the particles should belong to group B
groupA = hoomd.group.type(name='a-particles', type='A')
groupB = hoomd.group.type(name='b-particles', type='B')
NatomsA = len(groupA)
NatomsB = len(groupB)
print("Group A: ",NatomsA," particles (",100*NatomsA/Natoms,"%)")
print("Group B: ",NatomsB," particles (",100*NatomsB/Natoms,"%)")

################################################################
# 
# SET UP POTENTIAL
#
################################################################
print("\n\n\nSET UP POTENTIAL\n")
#Set up neighbor list
NeighborsListLJ = md.nlist.cell()

#Set Kob-Andersen potential
myLjPair=kobandersen.KApotential(NeighborsListLJ)
#These are two alternative ways of seeing the values we inserted
#print(myLjPair.pair_coeff.get_metadata())
#print(myLjPair.pair_coeff.values)



################################################################
# 
# SET UP ANALYZER
#
################################################################
print("\n\n\nSET UP ANALYZER\n")

#Name of the log
logname="test-output/observables-test_"+thermostat+str(nNVTsteps)+".txt"

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
# Initial FIRE minimization to avoid exploding the box,
# in case the initial configuration is very unstable
# 
################################################################

fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   hoomd.run(100)
print('FIRE minimization converged')

#This is what the Inherent Structure looks like
hoomd.dump.gsd(filename="./test-output/latticeIS.gsd", overwrite=True, period=None, group=hoomd.group.all(), time_step=0)


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
md.integrate.mode_standard(dt=dt)
md.update.zero_momentum(period=1,phase=0)


######## NVT ########
print(nNVTsteps," NVT steps with the ",thermostat," thermostat")
if nNVTsteps > 0 :

    if thermostat == 'NVT' :
       integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=0.5*TemperatureGoal, tau=tauT)
       print("Running at T=",TemperatureGoal*0.5)
       hoomd.run(nNVTsteps, quiet=False)
       integrator_nvt.set_params(tau=tauT, kT=TemperatureGoal)
       print("Running at T=",TemperatureGoal)
       hoomd.run(nNVTsteps, quiet=False)
        
    elif thermostat == 'MB' :
        stepsTauT = int(tauT/dt)
        integrator_nvt = md.integrate.nve(group=hoomd.group.all())
        for iterations in range(0,int(nNVTsteps/stepsTauT)):
            snap = system.take_snapshot()
            # each component is a gaussian of variance sigma, that's it.
            vel = np.random.normal(0,sqrt(TemperatureGoal), (Natoms,3))
            print('\nVelocities were rescaled by ',TemperatureGoal**0.5/np.std(vel),'\n')
            vel *= sqrt(TemperatureGoal)/np.std(vel,axis=0) #Riscalo le velocita` per essere sicuro che non ci sono effetti di taglia finita sulla varianza della distribuzione
            vel -=vel.mean(axis=0)
            snap.particles.velocity[:] = vel
            system.restore_snapshot(snap)
            hoomd.run(stepsTauT, quiet=False)
    else:
        print("In this tutorial the only implemented thermostats are NVT (Nose-Hoover) and MB (Anderson)\n")
        sys.exit()
    integrator_nvt.disable()
