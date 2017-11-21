### read-quench-record.py  ###
### reads a gsd file, for each frame, quench it to its local energy minimum using FIRE algorithm
### then record each frame to a gsd file (only the first frame contains all the information, other frames only contain the positions)
### input  : filename.gsd
### output : outName.gsd

import hoomd
from hoomd import md
import sys
import numpy as np
import module_filenamesHandler 
import gsd.pygsd
import gsd.hoomd
from subprocess import call 

import time

filename = sys.argv[1]
ftol = 1e-3
if len(sys.argv) > 2 :
    ftol = float(sys.argv[2])

NframesMax = 5000
pieceNumbe=0
if len(sys.argv) > 4 :
    NframesMax = int(sys.argv[3])
    pieceNumbe = int(sys.argv[4])
    
every = 1
Etol = 1e-5 ## not actually a constraint, becasue ftol is much more stringent.

## TODO: add "notes logaction in this 

#############################################################
rootname = module_filenamesHandler.get_rootname(filename)
suffix   = module_filenamesHandler.get_suffix(filename)

############################
### read the file (gsd) ####
with open(filename, 'rb') as flow:
    hoomdInputFlow = gsd.pygsd.GSDFile(flow)
    hoomdTraj = gsd.hoomd.HOOMDTrajectory(hoomdInputFlow);
    print 'using hoomdTraj.file.gsd_version=', hoomdTraj.file.gsd_version
    print 'using hoomdTraj.file.gsd_application=', hoomdTraj.file.application
    print 'using hoomdTraj.file.schema=', hoomdTraj.file.schema
    Nframes = hoomdTraj.file.nframes

    ## we assume these do not change over time:
    boxParams   = hoomdTraj.read_frame(0).configuration.box
    dim         = hoomdTraj.read_frame(0).configuration.dimensions
    Natoms      = hoomdTraj[0].particles.N
    atomIds     = np.arange(Natoms, dtype=int).transpose()
    atomTypes   = np.array(hoomdTraj[0].particles.typeid, dtype=int).transpose()

    print "there are Nframes=", Nframes
    print "we treat Nframes=", min(NframesMax, Nframes)
    tInit= (0+pieceNumbe)*NframesMax
    tFin = min((1+pieceNumbe)*NframesMax, Nframes)

    steps = np.zeros(Nframes)
    for t0 in range(tInit, tFin ,every):
        steps[t0] = hoomdTraj.read_frame(t0).configuration.step

outName = rootname+"_type=FIRE-traj_"+suffix[:-4]+"_piece="+str(pieceNumbe)+"_tInit="+str(tInit)+"_tFin="+str(tFin)+".gsd"
## erase the privous output file if it exists;
errMsg = call("rm "+outName, shell=True)


computeTime0 = time.time()
print "the number of frames to treat is ", Nframes , "  this is piece #",pieceNumbe, "   out of ", Nframes//NframesMax
for t0 in range(tInit, tFin ,every):
    step = steps[t0]
#    if t0<3 or t0>Nframes-every:
#        print "reading (firsts or lasts) steps: ", step

    notice_level=0
    hoomd.context.initialize('--notice-level='+str(notice_level))   
    # --mode=gpu 
    ## it is not absolutely necessary to use gpu.. 
    ## although it would help spped up things.
    system = hoomd.init.read_gsd(filename=filename, frame=t0)

    ################################################################################
    ########## Set up the interactions #############################################
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
    myLjPair.pair_coeff.set('A', 'A', epsilon=eps_AA, sigma=sig_AA, r_cut=r_cutoff*sig_AA)
    myLjPair.pair_coeff.set('A', 'B', epsilon=eps_AB, sigma=sig_AB, r_cut=r_cutoff*sig_AB)
    myLjPair.pair_coeff.set('B', 'B', epsilon=eps_BB, sigma=sig_BB, r_cut=r_cutoff*sig_BB)
    myLjPair.set_params(mode="shift")


    ### run FIRE minnimization/quench
    fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001, ftol=ftol, Etol=Etol)
    while not(fire.has_converged()):
       hoomd.run(500, quiet=True)
    print 'steps to converge the step #t0=',t0,' : ', hoomd.get_step() - step , ' FIRE time steps  and computeTime=',int(time.time() - computeTime0 )

    ### record quenched state.
    gsd_restart = hoomd.dump.gsd(filename=outName, group=hoomd.group.all(), period=None, overwrite=False, truncate=False)

#### FOR TESTING CONVERGENCE:  / DEBUG PURPOSES : 
#    pos = system.take_snapshot().particles.position
#    with open("StatesForCompare_"+str(step)+'_Etol='+str(Etol)+'_ftol='+str(ftol)+'_.dump', 'w') as outFlow:
#        np.savetxt(outFlow, pos)   
            


########### ftol = 1e-5:
#steps to converge:  10000.0
#steps to converge:  9000.0
#steps to converge:  15000.0
#steps to converge:  17500.0
#steps to converge:  13500.0
#steps to converge:  10000.0
#steps to converge:  11500.0
#steps to converge:  13500.0
#steps to converge:  9500.0
#steps to converge:  7500.0
#steps to converge:  12500.0
#steps to converge:  16000.0
#steps to converge:  12000.0
#steps to converge:  13500.0
#steps to converge:  19000.0
#steps to converge:  10500.0
#steps to converge:  13000.0
#steps to converge:  11500.0

########## ftol = 1e-3
#steps to converge:  6500.0
#steps to converge:  5500.0
#steps to converge:  6000.0
#steps to converge:  14500.0
#steps to converge:  6500.0
#steps to converge:  5500.0
#steps to converge:  8500.0
#steps to converge:  11000.0
#steps to converge:  7000.0
#steps to converge:  4000.0
#steps to converge:  9500.0
#steps to converge:  13500.0
#steps to converge:  10000.0
#steps to converge:  10000.0
#steps to converge:  14500.0
#steps to converge:  8000.0
#steps to converge:  10500.0


