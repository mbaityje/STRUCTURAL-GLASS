#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a trajectory from a .gsd file, and calculates the
# mean square displacement and self-intermediate scattering functions.
#
# To launch a simulation:
# python SingleSampleMeasurementsDynamic.py trajectory.gsd --dt=dt --every_forMemory=every_forMemory
#
# For example:
# python SingleSampleMeasurementsDynamic.py trajectory.gsd --dt=0.001 --every_forMemory=1
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import sys
from os import remove
import argparse
import numpy as np
import gsd.pygsd
import gsd.hoomd
import lib.module_measurements as med
from numba import jit

################################################################
#
# READ ARGUMENTS
# 
################################################################
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', #positional argument
                    nargs=1,
                    help='name of the .gsd trajectory we want to read'
)
parser.add_argument('--every_forMemory', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[1],
                    help='We only take one frame every every_forMemory frames'
)
parser.add_argument('--dt', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[0.0025],
                    help='dt of the MD dynamics (it is only used for the x axis in the figures and txt)'
)
parser.add_argument('-l','--label', #optional argument
                    nargs=1,
                    required=False,
                    default=[''],
                    help='label for distinguishing runs and continuations'
)
args = parser.parse_args()
filename=args.filename[0]
every_forMemory=args.every_forMemory[0]
dt=args.dt[0]
label=str(args.label[0])
height=0.36787944117144232159 #This is 1/e

#Wave vector for the self-intermediate scattring function
# k =[2 pi/L](n1,n2,n3) and permutations
n1=1
n2=3
n3=4


print("filename = ",filename)
print("dt = ",dt)
print("every_forMemory = ",every_forMemory)
print("label = ",label)
print("height for tau: ",height)
del parser

################################################################
#
# READ TRAJECTORY
# 
################################################################
with open(filename, 'rb') as flow:
    HoomdFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(HoomdFlow);
    s0=hoomdTraj.read_frame(0) #This is a snapshot of the initial configuration (frame zero)
    Natoms=s0.particles.N
    boxParams=s0.configuration.box
    L=np.float64(boxParams[0])
    time0=np.int64(s0.configuration.step)
    print("TIME_0:",time0)
    if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
        print('box dimensions are : ', boxParams[0])
        print('and are not supported (only isotropic systems supported)')
        raise SystemExit

    Nframes = len(hoomdTraj)
    trajDuration = Nframes/every_forMemory
    print('there are Nframes=', Nframes, 'in the file, but we only use trajDuration=',trajDuration, ' of them.')

    #Allocate memory
    trajectory = np.zeros((trajDuration,Natoms,3), dtype=np.float64)
    initialPositions=np.array(hoomdTraj[0].particles.position, dtype=np.float64)
    times=np.zeros(trajDuration, dtype=np.float64)
    msd=np.zeros(trajDuration, dtype=np.float64)
    Fk=np.zeros(trajDuration,dtype=np.float64)

    ################################################################
    #
    # SAVE TRAJECTORY AND CALCULATE MEAN SQUARE DISPLACEMENT
    #
    ################################################################
    for iframe in range(0, trajDuration, 1):
        ## we only look 1 frame every "every_forMemory" frame, to save memory: 
        trajectory[iframe] = np.array(hoomdTraj[iframe*every_forMemory].particles.position)
        time=np.int64(hoomdTraj[iframe*every_forMemory].configuration.step)-time0
        times[iframe]=time*every_forMemory*dt

        msd[iframe]=med.PeriodicSquareDistance(trajectory[iframe], initialPositions, L)/Natoms
        all_displacements=med.PeriodicDisplacement(trajectory[iframe], initialPositions, L)

        Fk[iframe]=med.ComputeFkt(n1, n2, n3, L, all_displacements)
    HoomdFlow.close()
    


################################################################
# 
# SAVE MSD and FkT
#
################################################################

#
# MSD and Fk(t)
#
#Remove preexisting output files to avoid conflicts
namemsd_txt='msd'+label+'.txt'
nameFkt_txt='Fkt'+label+'.txt'
namemsd_png='msd'+label+'.png'
nameFkt_png='Fkt'+label+'.png'
try:
    remove(namemsd_txt)
    remove(namemsd_png)
    remove(nameFkt_txt)
    remove(nameFkt_png)
except OSError:
    pass

#
# Text outputs
#
output_msd=np.column_stack((times, msd))
np.savetxt(namemsd_txt, output_msd,fmt='%g %.14g', header="#1)time step 2)msd")
output_Fk=np.column_stack((times, Fk))
np.savetxt(nameFkt_txt,output_Fk,fmt='%g %.14g', header="#n=("+str(n1)+","+str(n2)+","+str(n3)+")\n#1)time 2)Fk(t)")


#
# Figures
#
import matplotlib.pyplot as plt
plt.switch_backend('agg') #In order to be able to use pyplot without Xterm

#Figure of msd
plt.figure(1)
plt.loglog(times, msd)
plt.xlabel('t')
plt.ylabel('mean square displacement')
plt.title('Mean square displacement')
plt.grid(True)
plt.savefig(namemsd_png)

#Figure of Fk(t)
plt.figure(2)
axes = plt.gca()
axes.set_ylim([0,1])
line=np.empty(len(times)); line.fill(height)
plt.semilogx(times, Fk)
plt.semilogx(times, line, linewidth=0.1, color='black')
plt.xlabel('t')
plt.ylabel('self-intermediate scattering function')
plt.title('Self-intermediate scattering function')
plt.grid(True)
plt.savefig(nameFkt_png)



################################################################
# 
# CALCULATE TAU FROM THE Fk(t)
#
################################################################
tau=med.CalculateTau(Fk,0,Nframes-1, height=height,dt=dt,every_forMemory=every_forMemory)
print("tau=",tau)
nametau_txt='tau'+label+'.txt'
np.savetxt(nametau_txt, [tau],fmt='%.14g', header="#tau-- height= "+str(height))

