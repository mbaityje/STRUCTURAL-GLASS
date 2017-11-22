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
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)

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

def Bisect(t1, snap1, t2, snap2):
    #Due casi che non devono accadere
    assert(t2>=t1)
    if t2==t1:
        return;
    #Ora il caso in cui non est necessario comparare con e1 (t2-t1=1)
    e2=Minimize(snap2)
    if (t2-t1)==1:
        elist.append([t0+t2,e2])
        return
    #Ora comparo le due IS
    e1=Minimize(snap1)
    print("Le due eIS:",t1,e1," ,   ",t2,e2,"       diff=",e1-e2)
    #Se sono uguali salvo e finisco
    if(np.abs(e1-e2)<deltaE):
        elist.append([t0+t2,e2])
        return
    #Se sono differenti trovo il tempo intermedio (che est per forza
    #diverso per i check a inizio funzione) e biseziono sia il primo
    #che il secondo intervallo.
    else:
        t12=int(0.5*(t1+t2))
        snap12=system.take_snapshot()
        snap12.particles.position[:]=posizioni[t12]
        Bisect(t1,snap1,t12,snap12)
        Bisect(t12,snap12,t2,snap2)

################################################################
# 
# Bisection
#
################################################################
#List of energies
elist=[[t0,Minimize(snap_ini)]]
Bisect(0,snap_ini,Nframes-1,snap_final)



################################################################
# 
# Output
#
################################################################
if(ichunk==0):
    f_handle = file('elist.txt', 'w') #To write in overwrite mode
    np.savetxt(f_handle, elist,fmt='%d %.14g', header='#1)time 2)eIS')
else:
    f_handle = file('elist.txt', 'a') #To write in append mode
    np.savetxt(f_handle,elist,fmt='%d %.14g') 
f_handle.close()







