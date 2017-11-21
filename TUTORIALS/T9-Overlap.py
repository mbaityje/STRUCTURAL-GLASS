#!/usr/bin/python
################################################################
#
#
# DESCRIPTION
# - Start from a lattice configuration c1.
# - Reach a local inherent structure (IS) c2.
# - Warm up this inherent structure for t MD steps, reaching c3.
# - Reach IS c4.
# - Warm up for 100t MD steps, reaching c5.
# - Reach IS c6.
#
# Overlaps are calculated between all those configurations.
#
# To display help:
# python T9-Overlap.py --user="-h"
#
# To launch a simulation:
# python T9-Overlap.py --user="-t nNVTsteps"
#
# For example:
# python T9-Overlap.py --user="-t 10"
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
# FUNCTIONS THAT WOULD GO IN A SEPARATE MODULE
# 
################################################################

#Overlap with the configurations as input
def OverlapConfs(conf1, conf2,box_size):
   posizioni1=np.array(conf1.particles.position)
   posizioni2=np.array(conf2.particles.position)
   return OverlapPos(posizioni1,posizioni2,box_size)

#Overlap with the positions as input
def OverlapPos(posizioni1, posizioni2,box_size):
   dist=PeriodicDistance(posizioni1,posizioni2,box_size)
   return OverlapDist(dist,box_size)

#Overlap with the distance between confs as input
def OverlapDist(dist,box_size):
   delta=0.22 #half small particle diameter
   return np.where(dist<delta,1.,0.).sum()/Natoms


def PeriodicDistance(vec_a, vec_b, box_size):
# This function measures the distance between two lists of points, vec_a and vec_b,
# that can be vectors.
# box_size can be np.array([Lx,Ly,Lz]) or even just L
    delta = np.abs(vec_a - vec_b) #substraction and abs are done component by component
    delta = np.where(delta > 0.5 * box_size, delta - box_size, delta) #condition==True ? return second arg :otherwise return third
    return delta.sum(axis=-1)

 ################################################################
#
# SET UP THE SIMULATION CONTEXT
# 
################################################################
print("\n\n\nSET UP THE SIMULATION CONTEXT\n")
hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd



################################################################
#
# READ ARGUMENTS
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
                    default=[10],
                    help='number of NVT steps (default 0)'
)

args = parser.parse_args(more_arguments)
nNVTsteps=args.nNVTsteps[0]
print("nNVTsteps = ",nNVTsteps)


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



################################################################
# 
# COMPUTE CONFIGURATIONS
#
################################################################

#Initial Lattice c1
print("c1:")
c1 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=True, period=None, group=hoomd.group.all())

#Initial IS c2
print("c2:")
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   hoomd.run(100)
print('FIRE minimization converged')
c2 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=False, period=None, group=hoomd.group.all())

#First thermal configuration c3
print("c3:")
TemperatureGoal=1
tauT=1
dt=0.0025
md.integrate.mode_standard(dt=dt)
print(nNVTsteps," NVT steps")
integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), tau=tauT, kT=TemperatureGoal)
hoomd.run(nNVTsteps)
integrator_nvt.disable()
c3 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=False, period=None, group=hoomd.group.all())

#First thermal IS, c4
print("c4:")
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   hoomd.run(100)
print('FIRE minimization converged')
c4 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=False, period=None, group=hoomd.group.all())

#Second thermal configuration c5
print("c5:")
TemperatureGoal=1
tauT=1
dt=0.0025
md.integrate.mode_standard(dt=dt)
print(nNVTsteps," NVT steps")
integrator_nvt = md.integrate.nvt(group=hoomd.group.all(), tau=tauT, kT=TemperatureGoal)
hoomd.run(100*nNVTsteps)
integrator_nvt.disable()
c5 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=False, period=None, group=hoomd.group.all())

#Second thermal IS, c6
print("c6:")
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   hoomd.run(100)
print('FIRE minimization converged')
c6 = system.take_snapshot()
hoomd.dump.gsd(filename="./test-output/T9-Overlap.gsd", overwrite=False, period=None, group=hoomd.group.all())




################################################################
# 
# OVERLAPS
#
################################################################
box_size=np.array([system.box.Lx,system.box.Ly,system.box.Lz])



#Overlap with the configurations as input
print("q11 = ",OverlapConfs(c1,c1,box_size))
print("q12 = ",OverlapConfs(c1,c2,box_size))

#Overlap with the positions as input
pos1=np.array(c1.particles.position)
pos2=np.array(c2.particles.position)
pos3=np.array(c3.particles.position)
pos4=np.array(c4.particles.position)
pos5=np.array(c5.particles.position)
pos6=np.array(c6.particles.position)
print("q33 = ",OverlapPos(pos3,pos3,box_size))
print("q34 = ",OverlapPos(pos3,pos4,box_size))
print("q35 = ",OverlapPos(pos3,pos5,box_size))
print("Overlap between the two subsequent IS:")
print("q46 = ",OverlapPos(pos4,pos6,box_size))

#Overlap with the distances as input
dist56=PeriodicDistance(pos5,pos6,box_size)
print("Overlap between the last thermal configuration and its IS:", OverlapDist(dist56,box_size))


