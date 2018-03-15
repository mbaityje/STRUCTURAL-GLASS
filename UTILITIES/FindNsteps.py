#!/usr/bin/env python
###########################################################
#                                                         #
# This utility counts the number of frames in a gsd file. #
#                                                         #
###########################################################


from __future__ import print_function #for compatibility with python3.5
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import gsd.pygsd
import gsd.hoomd


filename = sys.argv[1] #Name of the file with the trajectory
iframe=0
if len(sys.argv)==3:
	iframe = int(sys.argv[2]) #Which frame
elif len(sys.argv)!=2:
    print("ERROR: Launch as:")
    print("python ",sys.argv[0]," configuration.gsd [ iframe=0]")
    sys.exit()

with open(filename, 'rb') as flow:
    HoomdFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
    s0=hoomdTraj.read_frame(iframe) #This is a snapshot of the initial configuration (frame zero)
    Natoms=s0.particles.N
    Nframes = len(hoomdTraj)

print(s0.configuration.step)
