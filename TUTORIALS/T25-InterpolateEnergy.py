#!/usr/bin/env python
################################################################
#
#
# DESCRIPTION
# This example reads two configurations, makes a linear interpolation of n points between them, 
# calculates the energy at each of those, and returns the energy profile along the interpolation.
# If both conf1 and conf2 are the same configuration, the program runs a few MD steps and compares
# initial and final configuration
# 
# 
# To display help:
# python T25-InterpolateEnergy.py --user="-h"
#
# To launch a simulation:
# python T25-InterpolateEnergy.py --user="conf1 conf2 --f1=frame1 --f2=frame2 -n=npoints"
#
# For example:
# python T25-InterpolateEnergy.py --user="./sample-states/rotenbergKA_T2.0_N1080.gsd ./sample-states/rotenbergKA_T2.0_N1080.gsd -n20"
# python T25-InterpolateEnergy.py --user="./sample-states/rotenbergKA_T2.0_N1080.gsd ./sample-states/rotenbergKA_T1.5_N1080.gsd -n20"
# 
################################################################

from __future__ import print_function 
import sys, argparse, copy
import numpy as np 
import hoomd, gsd.pygsd, gsd.hoomd
from hoomd import md
from lib import module_potentials as pot, module_measurements as med
np.set_printoptions(precision=15)

#Integrator parameters
alphaFIRE=0.99
ftolFIRE=1e-5
EtolFIRE=1e-10
wtolFIRE=1e-5
minstepsFIRE=100

hoomd.context.initialize() #Reads arguments from sys.argv and uses them automatically
hoomd.option.set_notice_level(0) #0: no notices - 10: lots of stuff
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename1', help='name of the first configuration file .gsd')
parser.add_argument('filename2', help='name of the second configuration file .gsd')
parser.add_argument('--f1', type=int, required=False, default=0, help='read the f1-th frame of filename1')
parser.add_argument('--f2', type=int, required=False, default=0, help='read the f2-th frame of filename2')
parser.add_argument('-n','--npoints', type=int, required=False, default=10, help='number of points in the interpolation')
args = parser.parse_args(more_arguments)

#Read First File
system=hoomd.init.read_gsd(args.filename1, frame=args.f1, time_step=0, restart=None)
snap=system.take_snapshot()
Natoms = snap.particles.N
L=snap.box.Lx
assert(snap.box.Lz==L and snap.box.Ly==L) #beware, these are floats, so there could be a false negative

#Read Second File
with open(args.filename2, 'rb') as flow:
	HoomdFlow = gsd.pygsd.GSDFile(flow)
	hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
	s2=hoomdTraj.read_frame(args.f2) 

# Check that the two configurations represent the same system
# Same number of particles
if s2.particles.N != Natoms:
	raise ValueError('The two configurations have different numbers of particles (%d and %d)'%(Natoms,s2.particles.N))
#Same particle ids
if all(snap.particles.typeid==s2.particles.typeid) == False:
	raise ValueError('The two configurations do not have the same particle indices')
#Same space dimensionality
if snap.box.dimensions != s2.configuration.dimensions:
	raise ValueError('The two configurations live in different spatial dimensions (%d and %d)'%(snap.box.dimensions, s2.configuration.dimensions))
#Same box size
if L==s2.configuration.box[0]==False:
	raise ValueError('Either the two boxes are different (Lx1=%g and Lx2=%g), or you should make a decent floating point comparison'%(snap.box.Lx, s2.configuration.box[0]))


#Neighborlist, Potential and Analyzer
NeighborsList = md.nlist.cell()
mypot="KAshort" if Natoms<500 else "KA"
potential=pot.LJ(NeighborsList,type=mypot)
analyzer = hoomd.analyze.log(filename=None, quantities=['temperature','potential_energy', 'kinetic_energy', 'momentum'], period=None, header_prefix = '#', overwrite=True, phase=0)

def Minimize(snap):
	system.restore_snapshot(snap)
	fire.cpp_integrator.reset()
	while not(fire.has_converged()):
		hoomd.run(100)
	eIS=analyzer.query('potential_energy')
	return eIS

#Initial and final positions
doMinimize=True
if doMinimize:
	fire=md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=alphaFIRE, ftol=ftolFIRE, Etol=EtolFIRE, wtol=wtolFIRE, min_steps=minstepsFIRE)
	integrator=md.integrate.nve(group=hoomd.group.all())
	snap.particles.velocity[:]=np.zeros((Natoms,3))
	Minimize(snap)
	eis=Minimize(snap)
	snap=system.take_snapshot()
	print("eis=",eis)
	pos1=np.array(snap.particles.position[:], dtype=np.float64)
	snap.particles.position[:]=s2.particles.position[:]
	system.restore_snapshot(snap)
	eis=Minimize(snap)
	print("eis=",eis)
	snap=system.take_snapshot()
	pos2=np.array(snap.particles.position[:], dtype=np.float64)
	integrator.disable()
else:
	pos1=np.array(snap.particles.position, dtype=np.float64)
	pos2=np.array(s2.particles.position, dtype=np.float64)


mode=md.integrate.mode_standard(dt=1e-20)
integrator=md.integrate.nve(group=hoomd.group.all())


#If the two configurations are the same, we evolve one of the two
if np.abs(pos1-pos2).sum()<1e-8:
	print('Running some MD steps...')
	hoomd.run(2)
	energy=analyzer.query('potential_energy')
	print('Energy1 = ',energy)
	mode.set_params(dt=0.001)
	hoomd.run(100)
	mode.set_params(dt=1e-20)
	energy=analyzer.query('potential_energy')
	print('Energy2 = ',energy)
	pos2=np.array( system.take_snapshot().particles.position[:],dtype=np.float64)

#Make interpolations
delta=1./args.npoints
alphas=np.arange(0, 1+delta, delta)
energies=[]
distances=[]
for alpha in alphas:
	interpPos=med.PeriodicInterpPoints(pos2, pos1, L, alpha)
	snap.particles.position[:]=interpPos
	system.restore_snapshot(snap)
	hoomd.run(2)
	energy=analyzer.query('potential_energy')
	print(energy)
	energies.append(energy)
	distances.append(med.PeriodicDistance(pos1,interpPos,L))

#Vector of all the particle displacements
disp=med.PeriodicDisplacement(pos1,pos2,L)

#
# FIGURES
#
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

#Barrier from the linear interpolation
fig, ax1 = plt.subplots()
ax1.set_xlabel('distance')
ax1.set_ylabel('energy')
ax1.plot(distances, energies)
ax1.tick_params(axis='x')
ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
ax2.set_xlabel('alpha')  # we already handled the x-label with ax1
ax2.plot(alphas, energies,'o',color='red',markersize=2)
ax2.tick_params(axis='x')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


#Vector plot of the Displacements between initial and final configuration
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
eps=1
ax.quiver(pos1[:,0], pos1[:,1], pos1[:,2], eps*disp[:,0], eps*disp[:,1], eps*disp[:,2])
ax.set_xlim([-L/2, L/2])
ax.set_ylim([-L/2, L/2])
ax.set_zlim([-L/2, L/2])
# ax.set_xlim([-L, L])
# ax.set_ylim([-L, L])
# ax.set_zlim([-L, L])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Displacement field between initial and final configuration')
plt.grid()
plt.show()




#Distribution of the displacements
sns.distplot(np.linalg.norm(disp,axis=1),bins=20)
plt.title('Distribution of displacements')
plt.show()





