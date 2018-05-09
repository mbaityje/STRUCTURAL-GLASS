#!/usr/bin/env python
################################################################
#
# This program calculates some correlation functions with high precision, and then fits them for sport.
# Then, it calculates the noise correlation function in a self-consistent way, starting from the
# calculated average correlations.
#
# DESCRIPTION
# This example:
# - reads a configuration from a .gsd
# - calculates forces
# - evolves the system and calculates force correlations at time t
# 
# The program calculates 4 autocorrelation functions:
# - velocity-velocity
# - force-force
# - force-velocity 
# - velocity-force (should be the same as the previous, with inverse sign)
# 
# 
# The relevant variables where one projects are
# P: (velocity of particle 1)
# F=dP/dt: (acceleration of particle 1)
# 
# Launch as:
# python T17-PreciseCorrelations.py --user="-N1000 -t10000 --dt=0.001 --temperature=1.5 --seed=0 --tauT=0.1"
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
parser.add_argument('-l','--label', required=False, default='', help='basename for the output files')
parser.add_argument('-s','--seed', type=int, required=False, default=-1, help='seed for random numbers. If negative we get it from /dev/urandom')
parser.add_argument('-t','--Ntraj', type=int, required=False, default=100, help='It is the total length of the run')
parser.add_argument('-n','--nchunk', type=int, required=False, default=10, help='Number of trajectories used to average the correlation')
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

Ncorr=int(args.Ntraj/args.nchunk)
print('Ncorr = ',Ncorr)

################################################################
#
# CREATE CONFIGURATION on an LxLxL lattice
#
################################################################
n=10  #Times the unit cell is repeated along each direction
a=1.259 #I choose this value so the rescaling is small (crashes with large rescale)
rho=0.5 #This is a very dilute liquid

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
system.box=hoomd.data.boxdim(volume=Natoms/rho)
print('density = ',rho)

################################################################
# FUNCTIONS
################################################################

def SaveMomentumForce(timestep):
	t=timestep-iniStep
	snap=system.take_snapshot()
	P[t]=np.array(snap.particles.velocity[:],dtype=np.float64)
	F[t]=np.array(snap.particles.acceleration[:],dtype=np.float64)
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

modeT=md.integrate.mode_standard(dt=args.dt)
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
print('Measure Correlations')
corrPP=np.zeros((args.nchunk,Ncorr))
corrFF=np.zeros((args.nchunk,Ncorr))
corrFP=np.zeros((args.nchunk,Ncorr))
corrPF=np.zeros((args.nchunk,Ncorr))
for ichunk in range(args.nchunk):
	t0=ichunk*Ncorr
	corrPP[ichunk]=np.array([np.mean([np.inner(P[t0][atom],P[time+t0][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
	corrFF[ichunk]=np.array([np.mean([np.inner(F[t0][atom],F[time+t0][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
	corrFP[ichunk]=np.array([np.mean([np.inner(F[t0][atom],P[time+t0][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
	corrPF[ichunk]=np.array([np.mean([np.inner(P[t0][atom],F[time+t0][atom]) for atom in range(Natoms)]) for time in range(Ncorr)])/3.
CPP   =np.mean(corrPP,axis=0,dtype=np.float64)
CFF   =np.mean(corrFF,axis=0,dtype=np.float64)
CFP   =np.mean(corrFP,axis=0,dtype=np.float64)
CPF   =np.mean(corrPF,axis=0,dtype=np.float64)
errCPP=stats.sem(corrPP,axis=0)
errCFF=stats.sem(corrFF,axis=0)
errCFP=stats.sem(corrFP,axis=0)
errCPF=stats.sem(corrPF,axis=0)





################################################################
# PLOT AND SAVE FIGURES
################################################################
xdata=np.arange(0,Ncorr)*args.dt

import matplotlib.pyplot as plt
plt.xlim((0, 1))
plt.errorbar(xdata, CPP/CPP[0], yerr=errCPP/CPP[0],errorevery=20, label='$\mathcal{C}^P_t/\mathcal{C}^P_0$', color='blue')
plt.errorbar(xdata, CFF/CFF[0], yerr=errCFF/CFF[0],errorevery=20, label='$\mathcal{C}^F_t/\mathcal{C}^F_0$', color='red')
plt.errorbar(xdata, 0*np.arange(0,Ncorr)          , yerr=None, label=None, color='black')
plt.xlabel('$t$')
plt.ylabel('$\mathcal{C}_t/\mathcal{C}_0$')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrPP-FF'+args.label+".png")
# plt.show()

plt.xlim((0, 1))
plt.ylim((-5.5, 5.5))
plt.errorbar(xdata, CPP/args.temperature, yerr=errCPP/args.temperature,errorevery=20, label='$\mathcal{C}^P_t/T$', color='blue')
# plt.errorbar(xdata, CFF/args.temperature, yerr=errCFF/args.temperature,errorevery=20, label='$\mathcal{C}^F_t/T$', color='red')
plt.errorbar(xdata, CFP/args.temperature, yerr=errCFP/args.temperature,errorevery=20, label='$\mathcal{C}^{FP}_t/T$', color='lawngreen')
plt.errorbar(xdata, CPF/args.temperature, yerr=errCPF/args.temperature,errorevery=20, label='$\mathcal{C}^{PF}_t/T$', color='lightgreen')
plt.errorbar(xdata, 0*np.arange(0,Ncorr)          , yerr=None, label=None, color='black')
plt.xlabel('$t$')
plt.ylabel('$\mathcal{C}_t/T$')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrOnT'+args.label+".png")
# plt.show()


################################################################
# FIT THOSE CORRELATION FUNCTIONS 
################################################################
from scipy.optimize import curve_fit

def f2(xdata, c1,lam1,c2,lam2):
	return c1*np.exp(-lam1*xdata)+c2*np.exp(-lam2*xdata)

xdata=np.arange(0,Ncorr)*args.dt
ydata=CPP
sigma=errCPP
my_fit=curve_fit(f2, xdata, CPP)
c1  =my_fit[0][0]
lam1=my_fit[0][1]
c2  =my_fit[0][2]
lam2=my_fit[0][3]

plt.errorbar(xdata, CPP, yerr=errCPP,errorevery=20, label='$\mathcal{C}^P_t/\mathcal{C}^P_0$', color='blue')
plt.errorbar(xdata, f2(xdata,c1,lam1,c2,lam2), label='fit $c_1e^{-\lambda_1t}+c_2e^{-\lambda_2t}$', color='lightblue')
plt.xlabel('$t$')
plt.ylabel('$\mathcal{C}_t/\mathcal{C}_0$')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrPP-fit'+args.label+".png")
# plt.show()


################################################################
# USE AVERAGED CORRELATORS TO CALCULATE THE NOISE FRICTION
################################################################
maxiter=100
Kold=np.copy(CFF)
Knew=np.zeros(Ncorr,dtype=np.float64)
xdata=np.arange(0,Ncorr)*args.dt

invT=np.float64(1./args.temperature)
for iter in range(maxiter):
	for n in range(Ncorr):
		Knew[n] = CFF[n] + invT * np.array([CFP[n-i]*Kold[i] for i in range(n)]).sum() * args.dt
	err=np.max(np.abs(Knew-Kold))
	# print(CFF[200],Kold[200], Knew[200])
	print("iter:",iter," err = ",err)
	if err<1e-10: break
	Kold[:]=Knew[:] #Curiosamente, se metto Kold=Knew, converge subito al risultato giusto

plt.plot(xdata, Knew,label='$\mathcal{K}_t$')
plt.plot(xdata, CFF,label='$\mathcal{C}^{F}_t$')
plt.xlabel('$t$')
plt.ylabel('Correlation')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrNoise-selfconsistent'+args.label+".png")
plt.show()


