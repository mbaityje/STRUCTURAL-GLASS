#!/usr/bin/env python
################################################################
#
# Same as T17, but for a Kob-Andersen mixture.
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
# The relevant variables where one projects are
# P: (velocity of particle 1)
# F=dP/dt: (acceleration of particle 1)
# 
# Launch as:
# python T18-OrthogonalCorrelationsKA.py --user="./sample-states/rotenbergKA_T10.0_N1080.gsd -N1080 -t10000 --dt=0.001 --temperature=10 --seed=0 --tauT=0.1 --extraThermalizing=True"
# python T18-OrthogonalCorrelationsKA.py --user="./sample-states/rotenbergKA_T10.0_N1080.gsd -N1080 -t80000 --dt=0.0025 --temperature=10 --seed=0 --tauT=0.1 --extraThermalizing=True --nchunk=200"
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
parser.add_argument('filename', help='name of the .gsd configuration we want to read')
parser.add_argument('-N','--Natoms', type=int, required=True, help='number of atoms in the system')
parser.add_argument('-l','--label', required=False, default='', help='basename for the output files')
parser.add_argument('-s','--seed', type=int, required=False, default=-1, help='seed for random numbers. If negative we get it from /dev/urandom')
parser.add_argument('-t','--Ntraj', type=int, required=False, default=100, help='It is the total length of the run')
parser.add_argument('-n','--nchunk', type=int, required=False, default=10, help='Number of trajectories used to average the correlation')
parser.add_argument('-T','--temperature', type=float, required=False, default=5, help='target Temperature')
parser.add_argument('--tauT', type=float, required=False, default=0.1, help='tau of the thermostat')
parser.add_argument('--dt', type=float, required=False, default=0.0025, help='dt for MD integration')
parser.add_argument('--extraThermalizing', type=bool, required=False, default=False, help='do some steps of thermalization before starting to measure')
args = parser.parse_args(more_arguments)
if args.seed>0: np.random.seed(args.seed)
print("filename = ",args.filename)
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
system = hoomd.init.read_gsd(filename=args.filename)
Natoms=len(system.particles)
assert(Natoms==args.Natoms)
print("Lattice with ",Natoms," particles")
rho=Natoms/system.box.get_volume()
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
potential=pot.LJ(NeighborsList,type="KA") if Natoms > 500 else pot.LJ(NeighborsList,type="KAshort")


modeT=md.integrate.mode_standard(dt=args.dt)
md.update.zero_momentum(phase=0, period=int(1./args.dt))

#Thermalize
if args.extraThermalizing:
	print('Thermalizing with NVT')
	integratorNVT = md.integrate.nvt(group=hoomd.group.all(), kT=args.temperature, tau=args.tauT)
	hoomd.run(int(10./args.dt), quiet=False)
	hoomd.dump.gsd(filename='./sample-states/rotenbergKA_T'+str(args.temperature)+'_N'+str(Natoms)+'.gsd', overwrite=True, truncate=True, period=None, time_step=0, group=hoomd.group.all())
	integratorNVT.disable()

#The trajectory
print('Measurement trajectory')
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
# plt.xlim((0, 1))
plt.errorbar(xdata, CPP/CPP[0], yerr=errCPP/CPP[0],errorevery=20, label='$\mathcal{C}^P_t/\mathcal{C}^P_0$', color='blue')
plt.errorbar(xdata, CFF/CFF[0], yerr=errCFF/CFF[0],errorevery=20, label='$\mathcal{C}^F_t/\mathcal{C}^F_0$', color='red')
plt.errorbar(xdata, 0*np.arange(0,Ncorr)          , yerr=None, label=None, color='black')
plt.xlabel('$t$')
plt.ylabel('$\mathcal{C}_t/\mathcal{C}_0$')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrPP-FF'+args.label+"_KA.png")
# plt.show()

# plt.xlim((0, 1))
# plt.ylim((-20, 20))
plt.errorbar(xdata, CPP/args.temperature, yerr=errCPP/args.temperature,errorevery=20, label='$\mathcal{C}^P_t/T$', color='blue')
# plt.errorbar(xdata, CFF/args.temperature, yerr=errCFF/args.temperature,errorevery=20, label='$\mathcal{C}^F_t/T$', color='red')
plt.errorbar(xdata, CFP/args.temperature, yerr=errCFP/args.temperature,errorevery=20, label='$\mathcal{C}^{FP}_t/T$', color='lawngreen')
plt.errorbar(xdata, CPF/args.temperature, yerr=errCPF/args.temperature,errorevery=20, label='$\mathcal{C}^{PF}_t/T$', color='lightgreen')
plt.errorbar(xdata, 0*np.arange(0,Ncorr)          , yerr=None, label=None, color='black')
plt.xlabel('$t$')
plt.ylabel('$\mathcal{C}_t/T$')
plt.grid(True)
plt.legend()
plt.savefig('./test-output/corrOnT'+args.label+"_KA.png")
# plt.show()




################################################################
# WITH TRUNCATED CORRELATION FUNCTIONS
################################################################
print('Measure Noise correlation functions')
maxiter=1000
Kold=np.copy(CFF)
Knew=np.zeros(Ncorr,dtype=np.float64)
xdata=np.arange(0,Ncorr)*args.dt

invT=np.float64(1./args.temperature)
for iter in range(maxiter):
	# Kold=np.where(np.abs(Kold)>nsigma*errCFP, Kold,0)
	for n in range(Ncorr):
		Knew[n] = CFF[n] + invT * np.array([CFP[n-i]*Kold[i] for i in range(n)]).sum() * args.dt
	err=np.max(np.abs(Knew-Kold))
	print("iter:",iter," err = ",err)
	if err<1e-10: break
	Kold[:]=Knew[:] #Curiosamente, se metto Kold=Knew, converge subito al risultato giusto

plt.plot(xdata, Knew,label='$\mathcal{K}_t$')
plt.plot(xdata, CFF,label='$\mathcal{C}^{F}_t$')
plt.xlabel('$t$')
plt.ylabel('Correlation')
plt.grid(True)
plt.legend()
plt.title('$\\rho = $%.1f,  $T = $%g,  $n = $%d,  $dt = $%g'%(rho,args.temperature,args.nchunk,args.dt))
plt.savefig('./test-output/corrNoise-selfconsistent'+args.label+"_KA.png")
plt.show()


