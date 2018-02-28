#!/usr/bin/python
################################################################
#
#
# DESCRIPTION
# This example creates a lattice configuration and melts it with a Nose-Hoover
# thermostat. The run stops when the system has thermalized, and saves the
# configuration.
#
# To display help:
# python  T3-CreateAndThermalize.py --user="-h"
#
# To launch a simulation:
# python  T3-CreateAndThermalize.py --user="-t nNVTsteps"
#
# For example:
# python T3-CreateAndThermalize.py --user="-t 5000"
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
# DEFINING FUNCTIONS
# 
################################################################

def IsAtEquilibrium34(array,verbose=False):
#Function takes observable in list, and compares its value in the
# third and fourth quarter of the list. Returns True if the two are compatible.
# Conditions for equilibration:
# -The two averages are statistically compatible
# -The two averages are within 1% relative difference
# -The statistical error of the fourth quarter is relatively <1%
# -The standard deviation of the fourth quarter is relatively <10% 
   quarter4=len(array)
   quarter3=(int)(3.*quarter4/4)
   quarter2=(int)(quarter4/2)
   part3=array[quarter2:quarter3]
   part4=array[quarter3:quarter4]
   mean3=part3.mean()
   mean4=part4.mean()
   std4=part4.std()
   err3=part3.std()/sqrt(quarter3-quarter2-1) #Statistical error of third quarter
   err4=std4/sqrt(quarter4-quarter3-1) #Statistical error of fourth quarter
   compatibility= abs(mean3-mean4)/sqrt(err3*err3+err4*err4) #Compatibility of the two measurements (in units of errorbars)
   relative_distance= abs(mean3-mean4)/mean4
   compatible= True if compatibility<1 else False #One std err apart
   relatively_close = True if relative_distance<0.01 else False #1% tolerance
   small_errorbar = True if err4/mean4<0.01 else False
   small_std = True if std4/mean4<0.1 else False
   if(verbose):
      print("quarter4: ",quarter4,"quarter3: ",quarter3,"quarter2: ",quarter2)
      print("Mean3: ", mean3," Err3: ",err3)
      print("Mean4: ", mean4," Err4: ",err4)
      print("Compatibility: ",compatibility," --> Compatible=",compatible)
      print("Relative Distance: ",relative_distance," --> Relatively close=",relatively_close)
      print("Relative Error Bar: ", err4/mean4," --> Small Error Bar=",small_errorbar)
      print("Relative Standard Deviation: ", std4/mean4," --> Small Standard Dev.=",small_std)
   return (compatible and relatively_close and small_errorbar and small_std) #Impose contraint both on statistical and on relative difference


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
# -t             number of NVT steps
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
                    default=[1000],
                    help='number of NVT steps (default 1000) between subsequent thermalization checks'
)
args = parser.parse_args(more_arguments)
nNVTsteps=args.nNVTsteps[0]
del parser
print(nNVTsteps)
assert(nNVTsteps>0)

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
uc = hoomd.lattice.unitcell(N = 5, #Number of particles in the unit cell: ___ Natoms=NxLxLxL ___
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
hoomd.dump.gsd(filename="./test-output/latticeThermalized.gsd", overwrite=True, period=None, group=hoomd.group.all(), time_step=0)
#Uncomment the next line if we want to save the trajectory
hoomd.dump.gsd(filename="./test-output/trajectoryThermalized.gsd", overwrite=True, period=1000, group=hoomd.group.all(), phase=0)

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
logname="test-output/observables-test_NVT"+str(nNVTsteps)+".txt"

#These are the observables we want to log
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']

print("In ",logname," we write ",analyzerManyVariables_quantities)

#Every how many integrations we want to log the observables
analyzer_period=100
analyzerManyVariables = hoomd.analyze.log(filename=logname, \
                                          quantities=analyzerManyVariables_quantities, period=analyzer_period, \
                                          header_prefix = '#This is a test output of the toy program'+sys.argv[0]+'\n#', \
                                          overwrite=True,
                                          phase=0)
analyzerManyVariables.disable() #I will enable it when the NVT part begins
 
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
hoomd.dump.gsd(filename="./test-output/latticeISthermalized.gsd", overwrite=True, period=None, group=hoomd.group.all(), time_step=0)


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
analyzerManyVariables.enable() #I will enable it when the NVT part begins


######## NVT ########
equilibrated=False
print(nNVTsteps," NVT steps with the Nose-Hoover thermostat at T=",TemperatureGoal)
integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)

while(False==equilibrated):
   print("Not yet at equilibrium, we run for ",nNVTsteps," more steps")
   hoomd.run(nNVTsteps, quiet=False)
   steplist,Tlist,Ulist,Klist=np.loadtxt(logname, skiprows=2, unpack=True, usecols=(0,1,3,4))
   Elist=np.array(Ulist)+np.array(Klist) #The list of the energy has the sum of potential and kinetic
   equilibrated= IsAtEquilibrium34(Tlist,verbose=False) and IsAtEquilibrium34(Elist,verbose=False)
print("Equilibrated!")

integrator_nvt.disable()


#This is what the Equilibrated Configuration looks like
hoomd.dump.gsd(filename="./test-output/equilibrated.gsd", overwrite=True, period=None, group=hoomd.group.all())
