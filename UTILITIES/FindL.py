#!/usr/bin/env python3
####################################################################
#                                                                  #
# This utility gives the size of the simulation box of a gsd file. #
#                                                                  #
####################################################################


from __future__ import print_function #for compatibility with python3.5
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import gsd.pygsd
import gsd.hoomd


if len(sys.argv)==2:
	filename = sys.argv[1] #Name of the file with the trajectory
else:
	print("ERROR: Launch as:")
	print("python ",sys.argv[0]," configuration.gsd")
	sys.exit()

with open(filename, 'rb') as flow:
	HoomdFlow = gsd.pygsd.GSDFile(flow)
	hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
	s0=hoomdTraj.read_frame(0) #This is a snapshot of the initial configuration (frame zero)
	Natoms=s0.particles.N

print("%.14g"%s0.configuration.box[0])

