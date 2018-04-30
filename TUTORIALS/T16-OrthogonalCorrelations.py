#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example:
# - reads a configuration from a .gsd
# - calculates forces
# - evolves the system and calculates force correlations at time t
# 
# The program calculates 3 autocorrelation functions:
# - velocity-velocity
# - force-force
# - noise-noise (i.e. the random force, i.e. the force in the orthogonal
# dynamics w.r. to the one defined through the projector operator)
# 
# 
# The relevant variables where one projects are
# P: (velocity of particle 1)
# F=dP/dt: (acceleration of particle 1)
# A: (force on particle 1) x k [k is a unit vector]
# B: (force on particle 1) x k
# 
# Launch as:
# python T16-OrthogonalCorrelations.py --user="-N1000 -t10000 --dt=0.001 --temperature=1.5 --seed=0 --tauT=0.1"
# 
################################################################




from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md #md is the subsection strictly about the properties of the system
import sys #this has some system tools like sys.exit() that can be useful
from os import remove #Used to remove the backup file at the end of the run
import argparse #for processing arguments
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #
from scipy import stats
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-N','--Natoms', type=int, required=True, help='number of atoms in the system')
parser.add_argument('-l','--label', required=False, default='thermalized', help='basename for the output files')
parser.add_argument('-s','--seed', type=int, required=False, default=-1, help='seed for random numbers. If negative we get it from /dev/urandom')
parser.add_argument('-t','--Ntraj', type=int, required=False, default=100, help='It is the total length of the run')
parser.add_argument('-T','--temperature', type=float, required=False, default=5, help='target Temperature')
parser.add_argument('--tauT', type=float, required=False, default=0.1, help='tau of the thermostat')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='dt for MD integration')
args = parser.parse_args(more_arguments)
if args.seed>0: np.random.seed(args.seed)
print("Natoms = ",args.Natoms)
print("seed = ",args.seed)
print("nNVTsteps = ",args.Ntraj)
print("T = ",args.temperature)
print("tauT = ",args.tauT)
print("dt = ",args.dt)
print("label = ",args.label)
beta=1./args.temperature
assert(args.Ntraj>0)
assert(args.temperature>0)
assert(args.tauT>0)
assert(args.dt>0 and args.dt<0.1)

Ncorr=int(args.Ntraj/10)
print('Ncorr = ',Ncorr)

################################################################
#
# CREATE CONFIGURATION on an LxLxL lattice
#
################################################################
n=10  #Times the unit cell is repeated along each direction
a=1.259 #I choose this value so the rescaling is small (crashes with large rescale)
rho=0.5 #This is a typical working density

try:
	system = hoomd.init.read_gsd(filename='./sample-states/rotenberg.gsd')
except RuntimeError:
	print("\n\nCREATE CONFIGURATION\n")
	#This is the unit cell we define
	#We need to define it with 5 elements, because the particle types are assigned at this stage
	#For a monodisperse mixture we could use the hoomd built-in unit cells
	uc = hoomd.lattice.unitcell(N = 1, #Number of particles in the unit cell
	                            a1 = [  a,   0,   0], #Basis vectors of the unit cell
	                            a2 = [  0,   a,   0],
	                            a3 = [  0,   0,   a],
	                            dimensions = 3,
	                            position = [[0, 0, 0]],
	                            type_name = ['A'],
	                            mass = [1.0],
	                            charge = [0.0],
	                            diameter = [1./2] #I believe that the diameter only plays a role in visualization
	);
	# lattice made of repeat of unit cell:
	system = hoomd.init.create_lattice(unitcell=uc, n=n)
Natoms=len(system.particles)
assert(Natoms==args.Natoms)
print("Lattice with ",Natoms," particles")
print("This is the",system.box)
print("Volume before rescaling: ",system.box.get_volume())
Box=hoomd.data.boxdim(volume=Natoms/rho)
system.box=Box
print("Volume after rescaling: ",system.box.get_volume())
print("This is the",system.box)
print('density = ',rho)

################################################################
# FUNCTIONS
################################################################

def SaveMomentumForce(timestep):
	t=timestep-iniStep
	snap=system.take_snapshot()
	P[t]=np.array(snap.particles.velocity[:])
	F[t]=np.array(snap.particles.acceleration[:])
	return



################################################################
# OBSERVABLES
################################################################
P=np.ndarray( (args.Ntraj, Natoms, 3), dtype=np.float64)
F=np.ndarray( (args.Ntraj, Natoms, 3), dtype=np.float64)

################################################################
# INTEGRATE
################################################################
NeighborsList = md.nlist.cell()
potential=pot.LJ(NeighborsList,type="LJmono")

modeT=md.integrate.mode_standard(dt=0.001)
md.update.zero_momentum(phase=0, period=int(1./args.dt))

#Thermalize
extraThermalizing=False
if extraThermalizing:
	print('Thermalizing with NVT')
	integratorNVT = md.integrate.nvt(group=hoomd.group.all(), kT=args.temperature, tau=args.tauT)
	hoomd.run(int(10./args.dt), quiet=False)
	hoomd.dump.gsd(filename='./sample-states/rotenberg.gsd', overwrite=True, truncate=True, period=None, time_step=0, group=hoomd.group.all())
	integratorNVT.disable()

#The trajectory
print('Measurement trajectory')
modeT.set_params(dt=args.dt)
integratorMeasure = md.integrate.nve(group=hoomd.group.all())
iniStep=hoomd.get_step()
# analyzerManyVariables = hoomd.analyze.log(filename=args.label+".txt", quantities=['temperature','potential_energy', 'kinetic_energy', 'momentum'], period=int(1./args.dt), header_prefix = '#', overwrite=True, phase=0)
# hoomd.dump.gsd(filename='./sample-states/trajectory'+args.label+'.gsd', overwrite=True, period=1, group=hoomd.group.all(),phase=-1)
callback=hoomd.analyze.callback(callback = SaveMomentumForce, period = 1, phase=0)
hoomd.run(args.Ntraj, quiet=False)
integratorMeasure.disable()



################################################################
# CALCULATE STANDARD CORRELATORS
################################################################
print('Velocity Correlation')
corrVel=np.array([np.mean([np.inner(P[0][atom],P[time][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
corrVelStd=np.array([stats.sem([np.inner(P[0][atom],P[time][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.

print('Force Correlation')
corrForce=np.array([np.mean([np.inner(F[0][atom],F[time][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
corrForceStd=np.array([stats.sem([np.inner(F[0][atom],F[time][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.





################################################################
# CALCULATE PROJECTED CORRELATOR
################################################################
def beta(my_n, my_Ntraj, my_F0, my_B, my_P0):
	mmax=my_Ntraj-my_n
	num = (my_F0[0:mmax]*my_B[n][0:mmax]).sum()
	den = np.square(my_P0[0:mmax]).sum()
	return num/den




print("Now calculate Projected Correlator")
klist=[ np.array([1,0,0], dtype=np.float64),
		np.array([0,1,0], dtype=np.float64),
		np.array([0,0,1], dtype=np.float64)]
taglist=range(1000)


nk=len(klist)
ntag=len(taglist)
for ik in range(nk):
	k=klist[ik]
	print('k = ',k)
	corrNoise=np.ndarray((ntag*nk,Ncorr),dtype=np.float64)
	for itag in taglist:
		print('itag = ',itag,'\t',end='')
		P0 = np.array([ P[t][itag].dot(k) for t in range(args.Ntraj) ])
		F0 = np.array([ F[t][itag].dot(k) for t in range(args.Ntraj) ])
		A0 = np.array([ F[t][itag].dot(k) for t in range(args.Ntraj) ])
		B=np.ndarray((Ncorr, args.Ntraj),dtype=np.float64)
		B[0] = np.array([ F[t][itag].dot(k) for t in range(args.Ntraj) ])

		print('B in Orthogonal Dynamics')
		for n in range(0, Ncorr-1):
			print('\rn=',n,end='')
			betan = beta(n, args.Ntraj, F0, B, P0) * args.dt
			mmax=args.Ntraj-n-1
			if 0:
				for m in range(mmax):
				# for m in range(args.Ntraj-1):
					B[n+1][m] = B[n][m+1] + betan*P0[m+1]
			else:
				# B[n+1][:-1] = B[n][1:]+betan*P0[1:]
				B[n+1][:mmax] = B[n][1:mmax+1]+betan*P0[1:mmax+1]



		print('\rCorrelation in Orthogonal Dynamics')
		for n in range(0, Ncorr):
			numel=args.Ntraj-n
			corrNoise[itag*nk+ik][n] = np.inner(A0[0:numel],B[n][0:numel])/numel


corrNoiseMean=np.zeros(Ncorr,dtype=np.float64)
corrNoiseStd =np.zeros(Ncorr,dtype=np.float64)


print('Output')
for n in range(0, Ncorr):
	corrNoiseMean[n] = np.mean(  corrNoise[:,n])
	corrNoiseStd [n] = stats.sem(corrNoise[:,n])
	print("%.4f"%(n*args.dt),corrNoiseMean[n],corrNoiseStd[n])





import matplotlib.pyplot as plt
plt.xlim((0, 1))
plt.ylim((-.25, 1))
plt.errorbar(np.arange(0,Ncorr)*args.dt, corrNoiseMean/corrNoiseMean[0], yerr=corrNoiseStd/corrNoiseMean[0], errorevery=20, label='Noise Mean')
plt.errorbar(np.arange(0,Ncorr)*args.dt, corrForce    /corrForce[0]    , yerr=corrForceStd/corrForce[0],     errorevery=20, label='Force')
plt.errorbar(np.arange(0,Ncorr)*args.dt, corrVel      /corrVel[0]      , yerr=corrVelStd  /corrVel[0],       errorevery=20, label='Velocity')
plt.errorbar(np.arange(0,Ncorr)*args.dt, 0*np.arange(0,Ncorr)          , yerr=None, label='Zero', color='black')
plt.xlabel('t')
plt.ylabel('corr')
plt.grid(True)
plt.legend()
# plt.savefig(filename+"_MSD.png")
plt.show()



