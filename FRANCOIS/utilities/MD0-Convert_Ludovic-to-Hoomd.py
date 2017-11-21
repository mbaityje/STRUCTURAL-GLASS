#!/usr/bin/env python2.7

################################################################################
# import hoomd and the md package
import hoomd
from hoomd import md
import sys
import numpy as np
import module_filenamesHandler

print sys.argv
headerLinesCount=0

print 'Usage: $ run convert_Ludovic-to-Hoomd.py  filename  special_tag   temperatureTimes100'

#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_1 9-1 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_2 9-2 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_3 9-3 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_4 9-4 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_5 9-5 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_6 9-6 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_7 9-7 41 &
#python  Convert_Ludovic-to-Hoomd.py DATA-Ludo/co_041_9_8 9-8 41 &


filename = sys.argv[1]
special_tag = sys.argv[2]
temperatureTimes100 = sys.argv[3]

#The proper Kob Anderson 80-20 (80% of A)
diameterA=1.0
diameterB=0.8
massA=1.0
massB=1.0
r_cutoff=2.5
eps_AA=1
eps_AB=1.5
eps_BB=0.5
sig_AA=1
sig_AB=0.8
sig_BB=0.88

################################################################################
######################### reading the files ####################################

print '\n\nopening ',filename,' assuming is it in the Ludovic\'s format'
Ep=[-42]
mysterious=[-42]
#xs,ys,zs = np.loadtxt(filename,skiprows=headerLinesCount, unpack=True)
xs,ys,zs,Ep,mysterious = np.loadtxt(filename,skiprows=headerLinesCount, unpack=True)

Natoms=len(xs)

atomTypes=np.zeros(Natoms)   
atomTypes[int(Natoms*0.8):]=1   ## the last 20% are type B (i.e. type 1)
#total of 4th column (probably E_p): -13815.8987
#which in Hoomd standards, div. by 2:  -6907.94935     (probably to not count twice)
#total of 5th column (???): 740.7077711


## reads a LAMMPS-style file
Lx=Ly=Lz=9.4103602888102831
extractedPositions = (np.array([xs-Lx/2.,ys-Ly/2.,zs-Lz/2.])).transpose()
################################################################################


################################################################################
##### creates an empty snapshot and then populate it with custom data. #########
hoomd.context.initialize('--notice-level=0') #'--mode=cpu') # '--mode=cpu')
myCustomSnapshot = hoomd.data.make_snapshot(Natoms, box=hoomd.data.boxdim(Lx=Lx,Ly=Ly,Lz=Lz), particle_types=['A','B'])
### hoomd.data.make_snapshot(Natoms, box, particle_types=['A'], bond_types=[], angle_types=[], dihedral_types=[], improper_types=[], dtype='float')
myCustomSnapshot.particles.position[:] = extractedPositions[:]
myCustomSnapshot.particles.typeid[:] = np.array(atomTypes[:], dtype=int)
diameters = np.ones(Natoms)
masses    = np.ones(Natoms)
diameters[myCustomSnapshot.particles.typeid[:]==0]=diameterA
diameters[myCustomSnapshot.particles.typeid[:]==1]=diameterB    ## autre tentativ: inverser A et B?
masses[   myCustomSnapshot.particles.typeid[:]==0]=massA
masses[   myCustomSnapshot.particles.typeid[:]==1]=massB
myCustomSnapshot.particles.diameter[:] = diameters[:]
myCustomSnapshot.particles.mass[:]     = masses[:]
print 'there is a percentage ', len(diameters[diameters==diameterA])*100.0/len(diameters), '%  of particles with a diameter ', diameterA
print 'there is a percentage ', len(diameters[diameters==diameterB])*100.0/len(diameters), '%  of particles with a diameter ', diameterB
### create the hoomd system from my custom-made snapshot:

system = hoomd.init.read_snapshot(myCustomSnapshot)

rootname= "KA_rho=12e-1_N="+str(Natoms)+"_T="+str(temperatureTimes100)+"e-2_tag="+special_tag
module_filenamesHandler.log_action(rootname, str(sys.argv)+'\n')
gsd_restart   = hoomd.dump.gsd(filename=rootname+"_type=restart-initial.gsd", truncate=True, period=None, group=hoomd.group.all(), phase=0)
gsd_restart.write_restart()

with open(rootname+"_type="+"notes.txt", 'a') as flow:
    flow.write("This set was created from \n"+filename+"\nIt was converted to Hoomd style gsd from a Ludovic Berthier format.\n")
    flow.write("From Ludo's file, we read a total E_p="+str(0.5*np.sum(Ep))+ '  in Hoomd notation\n')
    flow.write("From Ludo's file, we read a total E_k="+str(np.sum(mysterious))+ '  in Hoomd notation\n')
    flow.write("and we get, from standard analyzers (hoomd):  E_p=         E_k=        \n")


#gsd_restart2   = hoomd.dump.gsd(filename="restart-initial.gsd", truncate=True, period=None, group=hoomd.group.all(), phase=0)
#gsd_restart2.write_restart()

#gsd_restart.disable()
#gsd_restart2.disable()

############  THE END  #########################################################
################################################################################



print '\nThe job (converting files) is now done, but we run a bit of time steps to be sure we produced a valid restart file (and then we do NOT save that state) \n--------------------------------------------------------------------\n\n'
################################################################################
########## Set up the interactions #############################################

dt=0.005
r_cutoff=2.5
eps_AA=1 
eps_AB=1.5
eps_BB=0.5
sig_AA=1
sig_AB=0.8
sig_BB=0.88
## specify Lennard-Jones interactions between particle pairs
NeighborsListLJ = md.nlist.cell()
myLjPair = md.pair.lj(r_cut=r_cutoff, nlist=NeighborsListLJ)
myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA)
myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB)
myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB)
#myLjPair.set_params(mode="shift")

#### FIRE 
print '\nFIRE minimization starting \n'
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   print 'minimization not converged yet, runnning some more (or 1st time)'
   hoomd.run(5e3)
print '\nFIRE minimization converged\n'

snap_after_FIRE = system.take_snapshot()
#np.savetxt("snap_from-"+filename+"-after_FIRE.snap", snap_after_FIRE.particles.position)



## integrator parameters
md.integrate.mode_standard(dt=dt)


integrator_nve = md.integrate.nve(group=hoomd.group.all())
hoomd.run(1000)
integrator_nve.disable()

snap_after_NVE = system.take_snapshot()


