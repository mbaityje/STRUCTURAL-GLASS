#!/usr/bin/env python2.7

################################################################################
# import hoomd and the md package
import hoomd
from hoomd import md
import sys
import numpy as np

print sys.argv
headerLinesCount=9
## skip a few lines using skiprows=headerLinesCount  in np.loadtxt()

print 'Usage: $ run  XXXXXXXXX  filename  special_tag   temperatureTimes100'

foldername = sys.argv[1]
filename = sys.argv[2]
special_tag = sys.argv[3]
temperatureTimes100 = sys.argv[4]

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

## some other format where column 2 is missing (Tristan's monodipserse files)
#filename="fcc.1528620.dump"
#atomID, xs,ys,zs = np.loadtxt(filename,skiprows=headerLinesCount, unpack=True)
#Lx=8.0
#Ly=Lx
#Lz=Lx
#atomTypes=np.ones(Natoms)-1
#TemperatureGoal=0.58

atomID, atomTypes, xs,ys,zs = np.loadtxt(foldername+filename,skiprows=headerLinesCount, unpack=True)
atomTypes-=1    ## because they numbered them 1 and 2 instead of 0 and 1

Natoms=len(xs)
## reads a LAMMPS-style file
print '\n\nopening ',filename,' assuming is it in the LAMMPS format'
with open(foldername+filename, 'r') as f:
    lineNumber=0
    for line in f:
        lineNumber+=1
        if lineNumber<=headerLinesCount:
            if lineNumber in [1,3,5,9]:
                pass
            else:
                temp = np.fromstring(line, sep=' ')
                if lineNumber == 4 :
                    if Natoms != temp[0]: 
                        print 'the number of atoms in the header does not match the number of lines in file. \nExiting now\n\n'
                        raise SystemExit
                if lineNumber in [6,7,8]:
                    if temp[0] != 0 :
                        print 'the first coordinate of the box is not set to in x,y, or z. What to do?. \nExiting now\n\n'
                        raise SystemExit
                    if lineNumber==6:
                        Lx=temp[1]
                    if lineNumber==7:
                        Ly=temp[1]
                    if lineNumber==8:
                        Lz=temp[1]
            pass
        else:
            break
extractedPositions = (np.array([xs*Lx-Lx/2.,ys*Ly-Ly/2.,zs*Lz-Lz/2.])).transpose() 
################################################################################


################################################################################
##### creates an empty snapshot and then populate it with custom data. #########
hoomd.context.initialize() #'--mode=cpu') # '--mode=cpu')
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
    flow.write("This set was created from \n"+filename+"\nIt was converted to Hoomd style gsd from a LAMMPS dump file .\n")

gsd_restart2   = hoomd.dump.gsd(filename="restart.gsd", truncate=True, period=None, group=hoomd.group.all(), phase=0)
gsd_restart2.write_restart()

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
myLjPair.set_params(mode="shift")

#### FIRE 
print '\nFIRE minimization starting \n'
fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   print 'minimization not converged yet, runnning some more (or 1st time)'
   hoomd.run(5e3)
print '\nFIRE minimization converged\n'

snap_after_FIRE = system.take_snapshot()
#np.savetxt(foldername+"snap_from-"+filename+"-after_FIRE.snap", snap_after_FIRE.particles.position)



## integrator parameters
md.integrate.mode_standard(dt=dt)


integrator_nve = md.integrate.nve(group=hoomd.group.all())
hoomd.run(1000)
integrator_nve.disable()

snap_after_NVE = system.take_snapshot()


