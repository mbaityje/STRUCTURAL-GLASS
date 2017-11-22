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
parser.add_argument('filename', nargs=1, help='name of the initial state .gsd')
parser.add_argument('-d','--dt', nargs=1, type=float, required=False, default=[0.0025], help='dt for the MD dynamics')
parser.add_argument('--ichunk', nargs=1, type=int, required=True, help='index of the chunk')
parser.add_argument('--tchunk', nargs=1, type=int, required=True, help='number of MD steps in this chunk')
parser.add_argument('-l','--label', nargs=1, required=False, default=[''], help='label for distinguishing runs and continuations')
parser.add_argument('--temperature', nargs=1, type=float, required=True, help='In case we use an NVT thermostat')
parser.add_argument('--tauT', nargs=1, type=float, required=False, default=[0.1], help='tau of the thermostat')
args = parser.parse_args(more_arguments)

filename=args.filename[0]
dt=args.dt[0]
ichunk=args.ichunk[0]
tchunk=args.tchunk[0]
label=str(args.label[0])
temperature=args.temperature[0]
tauT=args.tauT[0]
del parser

print("filename: ",filename)
print("dt = ",dt)
print("ichunk = ",ichunk)
print("tchunk = ",tchunk)
print("temperature = ",temperature)
assert(temperature>=0)
print("tauT = ",tauT)

    


#Read configuration
restart_name_old    ="restartChunk"+str(ichunk-1)+label+".gsd"
restart_name_current="restartChunk"+str(ichunk)+label+".gsd"
if ichunk==0:
    system = hoomd.init.read_gsd(filename=filename,time_step=0)
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
hoomd.md.update.zero_momentum(period=int(1./dt),phase=0)


print("Now dynamics")
print("Current step:",hoomd.get_step())
print("Target step:", tchunk)
md.integrate.mode_standard(dt=dt)
integrator = md.integrate.nvt(group=hoomd.group.all(), kT=temperature, tau=tauT)
hoomd.md.update.zero_momentum(phase=-1)
hoomd.dump.gsd(filename="trajChunk"+str(ichunk)+label+".gsd", overwrite=True, period=1, group=hoomd.group.all(), phase=0)

#Now the MD steps
hoomd.run(tchunk)
hoomd.dump.gsd(filename=restart_name_current, group=hoomd.group.all(), truncate=True, period=None, phase=-1)
#try:
#    hoomd.run(tchunk)
#    hoomd.dump.gsd(filename=restart_name_current, group=hoomd.group.all(), truncate=True, phase=-1)
#except WalltimeLimitReached:
#    pass
integrator.disable()

