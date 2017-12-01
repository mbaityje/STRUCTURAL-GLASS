#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a trajectory and finds the related
# trajectory of the inherent structures.
#
################################################################

#MEMO:
#-Le posizioni devono essere in doppia precisione
#-Bisogna trovare il modo di fare accordare chunks successivi
#

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import argparse
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
import gsd.pygsd
import gsd.hoomd
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import module_measurements as med
import os.path
import copy

print("#++++++++++++++++#")
print("#----------------#")
print("# BisectChunk.Py #")
print("#----------------#")
print("#++++++++++++++++#")

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
parser.add_argument('filename', nargs=1, help='name of the trajectory file .gsd')
parser.add_argument('--ichunk', nargs=1, type=int, required=True, help='index of the chunk')
parser.add_argument('--tchunk', nargs=1, type=int, required=True, help='number of MD steps in this chunk')
parser.add_argument('-l','--label', nargs=1, required=False, default=[''], help='label for distinguishing runs and continuations')
parser.add_argument('--deltaE', nargs=1, type=float, required=False, default=[0.5], help='energy difference in order to consider different two configurations')
parser.add_argument('--skiprows', nargs=1, type=int, required=False, default=[0], help='how many rows we can skip when reading elist (we nee only the last one)')
args = parser.parse_args(more_arguments)

filename=args.filename[0]
ichunk=args.ichunk[0]
tchunk=args.tchunk[0]
deltaE=args.deltaE[0]
skiprows=args.skiprows[0]
label=str(args.label[0])
del parser

print("filename: ",filename)
print("ichunk = ",ichunk)
print("tchunk = ",tchunk)
print("deltaE = ",deltaE)


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
    print("Natoms = ",Natoms)
    boxParams=s0.configuration.box
    L=boxParams[0]
    if boxParams[0] != boxParams[1] or  boxParams[1]!= boxParams[2] or  boxParams[0]!= boxParams[2]:
        print('box dimensions are : ', boxParams[0])
        print('and are not supported (only isotropic systems supported)')
        raise SystemExit
    
    Nframes = len(hoomdTraj)
    print('There are Nframes=', Nframes, 'in the file.')
    assert(Nframes==tchunk)

    posizioni=[hoomdTraj[i].particles.position[:] for i in range(Nframes)]
    HoomdFlow.close()
        

################################################################
# 
# Initialize
#
################################################################
system = hoomd.init.read_gsd(filename=filename)
t0=hoomd.get_step()
print("Initial time:", t0)
snap_ini=system.take_snapshot()
snap_final=system.take_snapshot()
snap_final.particles.position[:]=posizioni[Nframes-1]

#If it's the first chunk, there is no list of energies.
#Otherwise, we open it and make sure that the time step is consistent.
if(ichunk>0):
    elist_old=np.loadtxt('elist.txt',skiprows=skiprows)
    assert(int(elist_old[len(elist_old)-1][0])==t0-1)

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
# Set analyzer
#
################################################################
print("\n\n\nSET UP ANALYZER\n")

analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
analyzer = hoomd.analyze.log(filename=None, quantities=analyzer_quantities, period=1)


################################################################
# 
# A couple of global variables:
#
################################################################
#Integrator
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.99, ftol=1e-2, Etol=1e-8, wtol=1e-2) #Do not reduce the precision
integrator = md.integrate.nve(group=hoomd.group.all())

################################################################
# 
# Function declaration
#
################################################################
def Minimize(snap):
    system.restore_snapshot(snap)
    fire.cpp_integrator.reset()
    while not(fire.has_converged()):
        hoomd.run(100)       
    eIS=analyzer.query('potential_energy')
    return eIS

def Bisect(t1, snap1, t2, snap2, e1=None, doTS=False):
    assert(t2>t1)
    # Energies of the inherent structures
    # If e1 was given as an argument, no need to caclulate it again
    if (e1==None):
        e1=Minimize(snap1)
    e2=Minimize(snap2)
    ediff=np.abs(e1-e2)
    print("Le due eIS:",t1,e1," ,   ",t2,e2,"       diff=",e1-e2)

    #If they are next to each other, I record them no matter their energy
    if (t2-t1)==1:
        if doTS and (ediff>deltaE): 
                eTS=NonLocalRidge(snap1, snap2, e1,e2)
        elist.append([t0+t2,e2])
        return

    #If they are the same I save and finish
    if(ediff<deltaE):
        elist.append([t0+t2,e2])
        return
    # If they are different, I find the intermediate time (that is necessarily
    # different because I made sure a few lines above), and I bisect both the first
    # and the second interval.
    else:
        t12=int(0.5*(t1+t2))
        snap12=system.take_snapshot()
        snap12.particles.position[:]=posizioni[t12]
        Bisect(t1,snap1,t12,snap12,e1, doTS=doTS)
        e12=Minimize(snap12)
        Bisect(t12,snap12,t2,snap2,e12, doTS=doTS)


def ConfBisect(snap1, snap2, e1, e2, L):
    #Maximum allowed distance between the configurations on the ridge
    dmax=0.002 #0.002 is about half the typical distance between confs at subsequent time steps w/ dt=0.0025
    dstart=0.1*dmax

    Natoms=snap_ini.particles.N
    pos1=np.array(snap1.particles.position, dtype=np.float64)
    pos2=np.array(snap2.particles.position, dtype=np.float64)
    dist12=med.PeriodicDistance(pos1,pos2,boxParams[0]).sum()/Natoms #the box is cubic

    count=0
    print("dist12=",dist12)
    while dist12>dstart:
        snap12=LinearConfInterpolation(snap1, snap2, L)
        e12=Minimize(snap12)
        print("e12=",e12)
        #If snap12 belongs to snap1, snap1=snap12
        if np.abs(e1-e12) <= 0.01: #Soglia scelta a cazzo, potrei mettere deltaE
            snap1=snap12
            e1=e12
            pos1=np.array(snap1.particles.position, dtype=np.float64)
        #If snap12 belongs to snap2, snap2=snap12
        elif np.abs(e2-e12) <= 0.01: #Soglia scelta a cazzo, potrei mettere deltaE
            snap2=snap12
            e2=e12
            pos2=np.array(snap2.particles.position, dtype=np.float64)
        #If snap12 does not belong to either, we throw a warning and change snap2
        else:
            print("NonLocalRidge: found an intermediate IS while searching the TS")
            snap2=snap12
            e2=e12
        #
        dist12=med.PeriodicDistance(pos1,pos2,boxParams[0]).sum()/Natoms
        print("dist12=",dist12)
        count+=1
        if count>10:
            print("NonLocalRidge ERROR: the interpolation bisection is not converging.")
    return snap1,snap2,snap12,e1,e2,e12,dist12


def NonLocalRidge(snap1, snap2, e1,e2):
    #Finds the Transition State through the algorithm proposed in Doliwa&Heuer, PRE 67 031506 (2003)
    print("Calculate TS now: |e1-e2|=",np.abs(e1-e2))

    snap1,snap2,snap12,e1,e2,e12,dist12=ConfBisect(snap1, snap2, e1, e2, np.float64(boxParams[0]))

    # ## Calculation of the gradient
    # grad=CalcolaGradiente(snap12)
    # g2=grad*grad
    # g2thres=0.01
    # niter=0
    # maxiter=10
    # while(g2>g2thres):
    #     snap1,snap2,dist12=ConfBisect()
    #     for iter in range(5):
    #         Minimizzo sia snap1 che snap2 per pochi passi
    #         dist12=med.PeriodicDistance()
    #         if dist12>dmax
    #             snap1,snap2,dist12=ConfBisect()
    #     grad=CalcolaGradiente()
    #     g2=grad*grad
    #     if niter>maxiter:
    #         print("IL CAZZO DI ALGORITMO DELLE SELLE DI MERDA NON CONVERGE MANCO SE LO PAGHI")
    #         Raise SystemExit
    #     niter++

    # ## Once the gradient is small, to steepest descent minimization of the square gradient
    # MinimizeSquareGradient()

    print("Interpolazione tra le due mi manda in una IS con e12=",Minimize(snap12))
    print("Mentre e1=",Minimize(snap1))
    print("Mentre e2=",Minimize(snap2))
    eTS=0
    return eTS

def LinearConfInterpolation(snap1, snap2, box_size):
#returns a linear interpolation between snap1 and snap2
#new_positions = (snap1 + snap2)/2
    pos1=np.array(snap1.particles.position,dtype=np.float64)
    pos2=np.array(snap2.particles.position,dtype=np.float64)
    snap12=system.take_snapshot()
    pos12=med.PeriodicIntermPoints(pos1,pos2,box_size)
    snap12.particles.position[:]=pos12
    return snap12

################################################################
# 
# Bisection
#
################################################################
#List of energies
elist=[[t0,Minimize(snap_ini)]]
Bisect(0,snap_ini,Nframes-1,snap_final, e1=None, doTS=True)

################################################################
# 
# Cleanup
#
################################################################
integrator.disable()

################################################################
# 
# Output
#
################################################################
if(ichunk==0):
    f_handle = file('elist.txt', 'w') #To write in overwrite mode
    np.savetxt(f_handle, elist,fmt='%d %.14g', header='1)time 2)eIS')
else:
    f_handle = file('elist.txt', 'a') #To write in append mode
    np.savetxt(f_handle,elist,fmt='%d %.14g') 
f_handle.close()







