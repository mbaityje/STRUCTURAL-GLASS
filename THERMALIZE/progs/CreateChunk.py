#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a configuration from a .gsd file. 
# The it runs the dynamics for a user-defined number of steps.
# Saves the thermal trajectory.
#
################################################################


from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import argparse
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import module_measurements as med
import os.path

################################################################
#
# READ ARGUMENTS
# 
################################################################
#Start hoomd
print("Initialize hoomd context\n")
simT=hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd



#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', help='name of the initial state .gsd')
parser.add_argument('-d','--dt', type=float, required=False, default=0.0025, help='dt for the MD dynamics')
parser.add_argument('--ichunk', type=int, required=True, help='index of the chunk')
parser.add_argument('--tchunk', type=int, required=True, help='number of MD steps in this chunk')
parser.add_argument('-l','--label', required=False, default='', help='label for distinguishing runs and continuations')
parser.add_argument('--temperature', type=float, required=False, help='In case we use an NVT thermostat')
parser.add_argument('--tauT', type=float, required=False, default=0.1, help='tau of the thermostat')
parser.add_argument('--thermostat', required=False, default='NVT', choices=['NVT','NVE'], help='What thermostat is used')
args = parser.parse_args(more_arguments)
del parser

print("filename: ",args.filename)
print("dt = ",args.dt)
print("ichunk = ",args.ichunk)
print("tchunk = ",args.tchunk)
print("thermostat = ",args.thermostat)
if args.thermostat=='NVT':
	assert(args.temperature>=0)
	print("T = ",args.temperature)
	print("tauT = ",args.tauT)


#Read configuration
restart_name_old    ="restartChunk"+str(args.ichunk-1)+args.label+".gsd"
restart_name_current="restartChunk"+str(args.ichunk)+args.label+".gsd"
if args.ichunk==0:
    system = hoomd.init.read_gsd(filename=args.filename,time_step=0)
else:
    system = hoomd.init.read_gsd(filename=restart_name_old)

Natoms = len(system.particles)
assert(system.box.Lx == system.box.Ly == system.box.Lz)
L=system.box.Lx
assert(system.box.get_volume()==L*L*L)

#More parameters
rho=Natoms/(system.box.get_volume())
print("Natoms = ",Natoms)
print("rho = ",rho)


################################################################
# 
# Set potential
#
################################################################
NeighborsListLJ = md.nlist.cell()
print(" *** KApotentialShort *** ")
myLjPair=pot.KApotentialShort(NeighborsListLJ)


################################################################
# 
# MD Dynamics
#
################################################################

#
#Write Initial configuration and make sure the center of mass is still
#
hoomd.md.update.zero_momentum(phase=-1)
hoomd.md.update.zero_momentum(period=int(1./args.dt),phase=0)


print("Now dynamics")
print("Current step:",hoomd.get_step())
print("Target step:", args.tchunk)
md.integrate.mode_standard(dt=args.dt)
if args.thermostat == 'NVT':
	integrator = md.integrate.nvt(group=hoomd.group.all(), kT=args.temperature, tau=args.tauT)
elif args.thermostat == 'NVE':
	integrator = md.integrate.nve(group=hoomd.group.all())
else:
	raise NotImplementedError('CreateChunk.py: The only implemented thermostats are NVE and NVT')

hoomd.md.update.zero_momentum(phase=-1)
hoomd.dump.gsd(filename="trajChunk"+str(args.ichunk)+args.label+".gsd", overwrite=True, period=1, group=hoomd.group.all(), phase=0)

#Now the MD steps
hoomd.run(args.tchunk)
hoomd.dump.gsd(filename=restart_name_current, group=hoomd.group.all(), truncate=True, period=None, phase=-1)
#try:
#    hoomd.run(args.tchunk)
#    hoomd.dump.gsd(filename=restart_name_current, group=hoomd.group.all(), truncate=True, phase=-1)
#except WalltimeLimitReached:
#    pass
integrator.disable()

