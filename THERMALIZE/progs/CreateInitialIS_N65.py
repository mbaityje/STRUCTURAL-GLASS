#!/usr/bin/python
################################################################
#
#
# DESCRIPTION
# This example creates a lattice configuration with N=65 and minimizes
# the energy. Save the initial configuration init.gsd, and the
# minimized one, initIS.gsd.
#
# To launch a simulation:
# python  CreateInitialIS.py
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import numpy as np #Handles some mathematical operations
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #
import sys 
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import  module_measurements as med #Here I have put the Kob-Andersen parameters


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


#Hard-coded variables:
Natoms=65
n_per_cell=Natoms
a=3.8 #I choose this value so the rescaling is small (crashes with large rescale)
rho=1.2 #This is a typical working density

print("Natoms = ",Natoms)
print("rho = ",rho)
print("a = ", a)

################################################################
#
# CREATE CONFIGURATION
#
################################################################

print("\n\n\nCREATE CONFIGURATION\n")

#We define one single unit cell of 65 elements.
#The unit cell is multiplied by a rescaling factor a (the size of the cell), and a security factor 0.9, because a is slightly larger than Lx,Ly,Lz at the chosen density. 
position=0.9*a*np.array([[-0.5  , -0.5  , -0.5  ],
       [-0.5  , -0.5  , -0.25 ],
       [-0.5  , -0.5  ,  0.   ],
       [-0.5  , -0.5  ,  0.25 ],
       [-0.5  , -0.25 , -0.5  ],
       [-0.5  , -0.25 , -0.25 ],
       [-0.5  , -0.25 ,  0.   ],
       [-0.5  , -0.25 ,  0.25 ],
       [-0.5  ,  0.   , -0.5  ],
       [-0.5  ,  0.   , -0.25 ],
       [-0.5  ,  0.   ,  0.   ],
       [-0.5  ,  0.   ,  0.25 ],
       [-0.5  ,  0.25 , -0.5  ],
       [-0.5  ,  0.25 , -0.25 ],
       [-0.5  ,  0.25 ,  0.   ],
       [-0.5  ,  0.25 ,  0.25 ],
       [-0.25 , -0.5  , -0.5  ],
       [-0.25 , -0.5  , -0.25 ],
       [-0.25 , -0.5  ,  0.   ],
       [-0.25 , -0.5  ,  0.25 ],
       [-0.25 , -0.25 , -0.5  ],
       [-0.25 , -0.25 , -0.25 ],
       [-0.25 , -0.25 ,  0.   ],
       [-0.25 , -0.25 ,  0.25 ],
       [-0.25 ,  0.   , -0.5  ],
       [-0.25 ,  0.   , -0.25 ],
       [-0.25 ,  0.   ,  0.   ],
       [-0.25 ,  0.   ,  0.25 ],
       [-0.25 ,  0.25 , -0.5  ],
       [-0.25 ,  0.25 , -0.25 ],
       [-0.25 ,  0.25 ,  0.   ],
       [-0.25 ,  0.25 ,  0.25 ],
       [ 0.   , -0.5  , -0.5  ],
       [ 0.   , -0.5  , -0.25 ],
       [ 0.   , -0.5  ,  0.   ],
       [ 0.   , -0.5  ,  0.25 ],
       [ 0.   , -0.25 , -0.5  ],
       [ 0.   , -0.25 , -0.25 ],
       [ 0.   , -0.25 ,  0.   ],
       [ 0.   , -0.25 ,  0.25 ],
       [ 0.   ,  0.   , -0.5  ],
       [ 0.   ,  0.   , -0.25 ],
       [ 0.   ,  0.   ,  0.   ],
       [ 0.   ,  0.   ,  0.25 ],
       [ 0.   ,  0.25 , -0.5  ],
       [ 0.   ,  0.25 , -0.25 ],
       [ 0.   ,  0.25 ,  0.   ],
       [ 0.   ,  0.25 ,  0.25 ],
       [ 0.25 , -0.5  , -0.5  ],
       [ 0.25 , -0.5  , -0.25 ],
       [ 0.25 , -0.5  ,  0.   ],
       [ 0.25 , -0.5  ,  0.25 ],
       [ 0.25 , -0.25 , -0.5  ],
       [ 0.25 , -0.25 , -0.25 ],
       [ 0.25 , -0.25 ,  0.   ],
       [ 0.25 , -0.25 ,  0.25 ],
       [ 0.25 ,  0.   , -0.5  ],
       [ 0.25 ,  0.   , -0.25 ],
       [ 0.25 ,  0.   ,  0.   ],
       [ 0.25 ,  0.   ,  0.25 ],
       [ 0.25 ,  0.25 , -0.5  ],
       [ 0.25 ,  0.25 , -0.25 ],
       [ 0.25 ,  0.25 ,  0.   ],
       [ 0.25 ,  0.25 ,  0.25 ],
       [-0.125, -0.125, -0.125]])
uc = hoomd.lattice.unitcell(N = n_per_cell,
                            a1 = [  a,   0,   0], #Basis vectors of the unit cell
                            a2 = [  0,   a,   0],
                            a3 = [  0,   0,   a],
                            dimensions = 3,
                            #64 particles in a grid, plus one at the center
                            position = position,
                            #An 80%-20% mixture
                            type_name = ['A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B',
                                         'A','A','A','A','B'],
                            charge = np.zeros(65),
                            mass = np.ones(65),
                            #I believe that the diameter only plays a role in visualization
                            diameter = [1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2,
                                        1./2, 1./2, 1./2, 1./2, 0.88/2
                                        ] 
);
# lattice made of repeat of unit cell:
system = hoomd.init.create_lattice(unitcell=uc, n=1)
assert(Natoms==len(system.particles))
print("Created a lattice with ",Natoms," particles")
print("This is the",system.box)
print("Volume before rescaling: ",system.box.get_volume())
Box=hoomd.data.boxdim(volume=Natoms/rho)
system.box=Box
print("Volume after rescaling: ",system.box.get_volume())
print("This is the",system.box)
hoomd.dump.gsd(filename="init.gsd", overwrite=True, period=None, group=hoomd.group.all())

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

#Set Kob-Andersen potential, adapted for small systems
print("KApotentialShort")
myLjPair=pot.KApotentialShort(NeighborsListLJ)


################################################################
# 
# Initial FIRE minimization to avoid exploding the box,
# in case the initial configuration is very unstable
# 
################################################################
print("FIRE minimization... ")
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.002)
integrator = md.integrate.nve(group=hoomd.group.all())
while not(fire.has_converged()):
   hoomd.run(100)
print('converged')
hoomd.dump.gsd(filename="initIS.gsd", overwrite=True, period=None, group=hoomd.group.all())


