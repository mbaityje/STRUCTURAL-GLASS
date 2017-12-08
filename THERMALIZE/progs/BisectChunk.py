#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example reads a trajectory and finds the related
# trajectory of the inherent structures.
#
# A search of the Transition States through Heuer's Non-Local
# Ridge Method is implemented, but I did not implement the final
# part of the algorithm (the minimization of the square gradient),
# so:
# ** run only with doTS=False **
#
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

potential=pot.LJ(NeighborsListLJ,type="KAshort") #myLJpair is now an attribute of potential. To call it: potential.GetLJpair()


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
# Function declaration
#
################################################################


#Returns system's energy if the analyzer is set
def PotEn():
    modeT.set_params(dt=1e-16)
    hoomd.run(2)
    U=analyzer.query('potential_energy')
    modeT.set_params(dt=0.0025)
    return U

def Minimize(snap):
    system.restore_snapshot(snap)
    fire.cpp_integrator.reset()
    while not(fire.has_converged()):
        hoomd.run(100)
    eIS=analyzer.query('potential_energy')
    return eIS

def MinimizeFewSteps(snap, nsteps):
    system.restore_snapshot(snap)
    fire.cpp_integrator.reset()
    hoomd.run(nsteps)
    snapnew=system.take_snapshot()
    return snapnew

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


def ConfBisect(snap1, snap2, eis1, eis2, L, dmax=0.002):
    #0.002 is about half the typical distance between confs at subsequent time steps w/ dt=0.0025

    assert(np.abs(eis1-eis2)>deltaE)
    dstart=0.1*dmax
    Natoms=snap_ini.particles.N
    pos1=np.array(snap1.particles.position, dtype=np.float64)
    pos2=np.array(snap2.particles.position, dtype=np.float64)
    dist12=med.PeriodicDistance(pos1,pos2,boxParams[0]).sum()/Natoms #the box is cubic
    print("dist12=",dist12)
    snap12=LinearConfInterpolation(snap1, snap2, L)
    pos12=np.array(snap12.particles.position, dtype=np.float64)
    eis12=Minimize(snap12)

    print(" eis12=",eis12)
    count=0
    maxcount=10
    while dist12>dstart:
        #If snap12 belongs to snap1, snap1=snap12
        if np.abs(eis1-eis12) <= deltaE:
            snap1=snap12
            eis1=eis12
            pos1=np.array(snap1.particles.position, dtype=np.float64)
        #If snap12 belongs to snap2, snap2=snap12
        elif np.abs(eis2-eis12) <= deltaE:
            snap2=snap12
            eis2=eis12
            pos2=np.array(snap2.particles.position, dtype=np.float64)
        #If snap12 does not belong to either, we throw a warning and change snap2
        else:
            print("NonLocalRidge: found an intermediate IS while searching the TS")
            snap2=snap12
            eis2=eis12
        #
        dist12=med.PeriodicDistance(pos1,pos2,boxParams[0]).sum()/Natoms
        print("dist12=",dist12)
        count+=1
        assert(np.abs(eis1-eis2)>deltaE)
        snap12=LinearConfInterpolation(snap1, snap2, L)
        eis12=Minimize(snap12)
        if count>maxcount:
            print("NonLocalRidge ERROR: the interpolation bisection is not converging.")
    return snap1,snap2,snap12,eis1,eis2,eis12,dist12


def NonLocalRidge(snap1, snap2, eis1,eis2):
    """Finds the Transition State through the algorithm proposed in Doliwa&Heuer, PRE 67 031506 (2003)"""
    print("Calculate TS now: |e1-e2|=",np.abs(eis1-eis2))
    dmax=0.002 #0.004 is about the typical distance between confs at subsequent time steps w/ dt=0.0025

    snap1,snap2,snap12,eis1,eis2,eis12,dist12=ConfBisect(snap1, snap2, eis1, eis2, np.float64(boxParams[0]),dmax=dmax)



    ## Calculation of the gradient

    e1,G1=potential.CalculateGradient(snap1)
    e2,G2=potential.CalculateGradient(snap2)
    Gsq1=np.square(G1).sum()
    Gsq2=np.square(G2).sum()
    Gsq=0.5*(Gsq1+Gsq2) #Questa variabile ha senso solo dopo ConfBisect, che garantisce che le configurazioni siano vicine
    eis1=Minimize(snap1)
    eis2=Minimize(snap2)

    GsqThres=0.02
    niter=0
    maxiter=200
    print("Nonlocal Ridge, caratteristiche configurazione termica:")
    print("e1=",e1)
    print("e2=",e2)
    print("Gsq1=",Gsq1)
    print("Gsq2=",Gsq2)
    print("Struttura inerente:")
    print("eis1=",eis1/Natoms)
    print("eis2=",eis2/Natoms)
    nsteps=25
    while(Gsq>GsqThres):

        snap1=MinimizeFewSteps(snap1, nsteps)
        snap2=MinimizeFewSteps(snap2, nsteps)
        pos1=np.array(snap1.particles.position,dtype=np.float64)
        pos2=np.array(snap2.particles.position,dtype=np.float64)
        dist12=med.PeriodicDistance(pos1,pos2,boxParams[0]).sum()/Natoms #the box is cubic
        # if dist12>dmax:
        snap1,snap2,snap12,eis1,eis2,eis12,dist12=ConfBisect(snap1, snap2, eis1, eis2, np.float64(boxParams[0]), dmax=dmax)

        e1,G1=potential.CalculateGradient(snap1)
        e2,G2=potential.CalculateGradient(snap2)
        e12,G12=potential.CalculateGradient(snap12)
        # eis1=Minimize(snap1)
        # eis2=Minimize(snap2)
        Gsq1=np.square(G1).sum()
        Gsq2=np.square(G2).sum()
        Gsq12=np.square(G12).sum()
        Gsq=0.5*(Gsq1+Gsq2)
        print("e1=",e1,"\te2=",e2)
        print("eis1=",eis1/Natoms,"\teis2=",eis2/Natoms)
        print("Gsq1=",Gsq1,"  Gsq2=",Gsq2,"  Gsq12=",Gsq12,"  Gsq=",Gsq,"\n")
        if np.abs(eis1-eis2)<deltaE:
            print("The two configurations are falling into the same IS. Abort.")
            raise SystemExit

        if niter>maxiter:
            print("IL CAZZO DI ALGORITMO DELLE SELLE DI MERDA NON CONVERGE MANCO SE LO PAGHI")
            raise SystemExit
        niter+=1

    snap1,snap2,snap12,eis1,eis2,eis12,dist12=ConfBisect(snap1, snap2, eis1, eis2, np.float64(boxParams[0]),dmax=dmax)
    eTS=e12
    e12,G12=potential.CalculateGradient(snap12)
    Gsq12=np.square(G12).sum()
    print("Ora bisogna fare una minimizzazione del gradiente quadro, a partire da snap12")
    print("Per la cronaca, l'energia di snap12 e` e12=",e12)
    print("Quella della struttura inerente e` eis12=",eis12)
    print("Il gradiente della configurazione termica e` Gsq12=",Gsq12)
    raise SystemExit
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

def CalcF():
    #Calculates the forces of the system through the integrator
    #Mind that if the dynamics were not run, this gives zero
    F=np.array([system.particles[i].net_force for i in range(Natoms)])
    return F/Natoms

def CalcAcc(snapshot):
    #Returns the accelerations of the snapshot
    #Mind that some snapshots don't have information on the accelerations
    acc=np.array(snapshot.particles.acceleration,dtype=np.float64)
    return acc/Natoms
  
################################################################
# 
# Bisection
#
################################################################
fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.99, ftol=1e-5, Etol=1e-10, wtol=1e-5) #Do not reduce the precision
integrator = md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)

print("Proprieta` dello snapshot")
snap=system.take_snapshot()
F=CalcF()
F2=np.square(F).sum()
A=CalcAcc(snap)
A2=np.square(A).sum()
e,G=potential.CalculateGradient(snap)
G2=np.square(G).sum()

print("F2=",F2)
print("A2=",A2)
print("G2=",G2)

print("Proprieta` di snap_ini (per il quale non ho runnato dinamica)")
system.restore_snapshot(snap_ini)
F=CalcF()
F2=np.square(F).sum()
A=CalcAcc(snap_ini)
A2=np.square(A).sum()
e,G=potential.CalculateGradient(snap_ini)
G2=np.square(G).sum()

print("F2=",F2)
print("A2=",A2)
print("G2=",G2)


#List of energies
elist=[[t0,Minimize(snap_ini)]]
Bisect(0,snap_ini,Nframes-1,snap_final, e1=None, doTS=False) #The search of the transition state (TS=True) does not work.


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







