## this code creates a kob-andersen mixture (not sure all potentials are correctly set)
## then it quenches the state, 
## then it runs it for a bit of time
## then it records it as a complete restart file (gsd, hoomd natural)

import hoomd
from hoomd import md
import numpy as np

# initialize the execution context
hoomd.context.initialize('--mode=cpu')

# unit cell size in our units (controls density)
#phi=0.77
#phi = (4*Va+1*Vb)/a**3
#a=((4*Va+1*Vb)/phi)**(1/3)
L_latticeSize=4
a=2.0

######################################################################################################
############### init for unit cell, lattice, atom types, init positions, etc. ########################
# initial positions of the 5 atoms in the unit cell 
initPos = a * np.array([[0,0,0], [0.5, 0.5, 0.5], [0.2,0.2,0.2], [0.7,0.7,0.7], [0.2,0.5,0.7]])
# custom-made unit cell : 
mySimpleCubicLattice = hoomd.lattice.unitcell(N = 5,
    a1 = [a,0,0],    a2 = [0,a,0],    a3 = [0,0,a],     ## shape of the unit cell
    dimensions = 3,
    position = initPos,
    type_name = ['A1','A2','A3', 'A4','C'],     # unfortunately it seems that Identical atoms are not allowed.. maybe I can create first all of type A, then modify type A->C ?
    diameter = [1.0, 1.0, 1.0, 1.0, 1.5]
#                            moment_inertia = [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]],
#                            orientation = [[0.707, 0, 0, 0.707], [1.0, 0, 0, 0]]
);

# lattice made of repeat of unit cell:
system = hoomd.init.create_lattice(unitcell=mySimpleCubicLattice, n=L_latticeSize)

# specify Lennard-Jones interactions between particle pairs
myNeighborsList = md.nlist.cell()
myLjPair = md.pair.lj(r_cut=3.0, nlist=myNeighborsList)
## what about the  singularity at the cutoff ??
eps_AA=1.0
sig_AA=1.0
eps_AB=1.5
sig_AB=0.8
myLjPair.pair_coeff.set('A1', 'A1', epsilon=eps_AA, sigma=sig_AA, r_on=2.4)  # use r_on < r_cut so that the xplor option also smoothes the Froce, not just the potential V
myLjPair.pair_coeff.set('A1', 'A2', epsilon=eps_AA, sigma=sig_AA, r_on=2.4)
myLjPair.pair_coeff.set('A1', 'A3', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A1', 'A4', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A2', 'A3', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A2', 'A4', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A3', 'A4', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A1', 'C' , epsilon=eps_AB, sigma=sig_AB, r_on=2.4)
myLjPair.pair_coeff.set('A2', 'A2', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A2', 'C' , epsilon=eps_AB, sigma=sig_AB, r_on=2.4)
myLjPair.pair_coeff.set('A3', 'A3', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A3', 'C' , epsilon=eps_AB, sigma=sig_AB, r_on=2.4)
myLjPair.pair_coeff.set('A4', 'A4', epsilon=eps_AA, sigma=sig_AA, r_on=2.4) 
myLjPair.pair_coeff.set('A4', 'C' , epsilon=eps_AB, sigma=sig_AB, r_on=2.4)
myLjPair.pair_coeff.set('C' , 'C' , epsilon=0.5, sigma=0.88     , r_on=2.4)

#myLjPair.set_params(mode="shift") # A constant shift is applied to the entire potential so that it is 0 at the cutoff
myLjPair.set_params(mode="xplor") #  A smoothing function is applied to gradually decrease both the force and potential to 0 at the cutoff when ron < rcut, and shifts the potential to 0 at the cutoff when ron >= rcut.


######################################################################################################
############### initial FIRE minimization to avoid exploding the box . ###############################

fire=hoomd.md.integrate.mode_minimize_fire(group=hoomd.group.all(), dt=0.001) # , ftol=1e-2, Etol=1e-7)
while not(fire.has_converged()):
   hoomd.run(10000)
print 'FIRE minimization converged'


######################################################################################################
############### running to thermalize the system a bit. ##############################################

# integrate with brownian dynamics
md.integrate.mode_standard(dt=0.001)
myIntegrator = md.integrate.brownian(group=hoomd.group.all(), kT=0.001, seed=987)

## time for the auto-tuner to find its values 
hoomd.run(1000)     

dcd_trajectory      = hoomd.dump.dcd(filename="trajectory_KobAndersen.dcd", period=1e3, group=hoomd.group.all(), phase=0,overwrite=False)
hoomd.dump.gsd(filename="init_KobAndersen.gsd", truncate=True, period=1e4, group=hoomd.group.all(), phase=0)   ## restart-like dump
##running for some time to thermalize the intiial condition, which is quenched at T=0 initially
hoomd.run(1e4)  

my_snap = system.take_snapshot()
np.savetxt("final_snapshot_positions.snap", my_snap.particles.position)
#np.savetxt("final_snapshot_positions.snap", my_snap.particles)





