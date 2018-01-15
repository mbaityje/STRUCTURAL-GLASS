#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This program reads a configuration from a .gsd file. 
# The it runs the dynamics for a user-defined number of steps.
# Saves e(t), eIS(t), q(0,t), qIS(0,t)
#
################################################################

from __future__ import print_function #for compatibility with python3.5
import hoomd #hoomd is the MD package from the Glotzer group
from hoomd import md
import argparse
import sys #this has some system tools like sys.exit() that can be useful
import numpy as np #Handles some mathematical operations
from lib import module_potentials as pot #Here I have put the Kob-Andersen parameters
from lib import module_measurements as med
import os.path

################################################################
#
# READ ARGUMENTS
# 
################################################################
#Start hoomd
print("Initialize hoomd context\n")
simT=hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user()

#Read parameters
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', nargs=1, help='name of the initial state .gsd')
parser.add_argument('-d','--dt', nargs=1, type=float, required=False, default=[0.0025], help='dt for the MD dynamics')
parser.add_argument('-l','--label', nargs=1, required=False, default=[''], help='label for distinguishing runs and continuations')
parser.add_argument('--temperature', nargs=1, type=float, required=True, help='In case we use an NVT thermostat')
parser.add_argument('--nsteps', nargs=1, type=int, required=True, help='Total number of steps')
parser.add_argument('--ntbar', nargs=1, type=int, required=False, default=[10], help='Number of times')
parser.add_argument('--seed', nargs=1, type=int, required=False, default=[-1], help='Seed for random numbers (default: 12345)')
parser.add_argument('--iframe', nargs=1, type=int, required=False, default=[0], help='Frame to read from the gsd file (default is 0)')
parser.add_argument('--tauT', nargs=1, type=float, required=False, default=[0.1], help='tau of the thermostat')
parser.add_argument('--irep', nargs=1, type=int, required=False, default=[0], help='replica index, for running many times the same trajectory')
args = parser.parse_args(more_arguments)

filename=args.filename[0]
dt=args.dt[0]
label=str(args.label[0])
temperature=args.temperature[0]
tauT=args.tauT[0]
nsteps=args.nsteps[0]
iframe=args.iframe[0]
ntbar=args.ntbar[0]
irep=args.irep[0]
seed=args.seed[0]

print("filename: ",filename)
print("dt = ",dt)
assert(dt>0 and dt<0.01)
print("nsteps = ",nsteps)
assert(nsteps>=0)
print("temperature = ",temperature)
assert(temperature>=0)
print("tauT = ",tauT)
assert(tauT>=0)
    

# Read configuration
system = hoomd.init.read_gsd(filename=filename,time_step=0, frame=iframe)
Natoms = len(system.particles)
c0=system.take_snapshot()
# Initialize c0 with custom velocities
if seed>=0:
	np.random.seed(seed)
	vel = np.random.normal(0,np.sqrt(temperature), (Natoms,3))
	vel *= np.sqrt(temperature)/np.std(vel,axis=0)
	vel -=vel.mean(axis=0)
	c0.particles.velocity[:] = vel


assert(system.box.Lx == system.box.Ly == system.box.Lz)
L=system.box.Lx
assert(system.box.get_volume()==L*L*L)


################################################################
# 
# Set potential
#
################################################################
NeighborsListLJ = md.nlist.cell()
potential=pot.LJ(NeighborsListLJ,type="KAshort")

################################################################
# 
# Set analyzer
#
################################################################
analyzer_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']
analyzer = hoomd.analyze.log(filename=None, quantities=analyzer_quantities, period=1)



################################################################
# 
# MD Dynamics
#
################################################################

#List of times
def ListaLogaritmica(x0,xn,n,ints=False,addzero=False):
    assert(xn>x0)
    assert(x0>0)
    n=np.int64(n)
    y0=np.log(x0)
    yn=np.log(xn)
    delta=np.float64(yn-y0)/(n-1)
    listax=np.exp([y0+i*delta for i in range(n)])
    if ints:
        listax=np.unique(np.round(listax,0).astype(int))
    if addzero:
        listax=np.insert(listax,0,0)
    return listax

tbar=ListaLogaritmica(1,nsteps,ntbar,ints=True,addzero=True)
print("Tiempos:",tbar)
EISlist=[]
Elist=[]
qlist=[1]
qISlist=[1]

#Initial Thermal Energy
system.restore_snapshot(c0)
md.integrate.mode_standard(dt=1e-16)
integrator_nvt = md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)
E0=analyzer.query('potential_energy')
# print("E0=",E0, )
integrator_nvt.disable()
Elist.append(E0)

# Initial IS
fire=hoomd.md.integrate.mode_minimize_fire(dt=dt, alpha_start=0.99, ftol=1e-5, Etol=1e-10, wtol=1e-5)
integrator_nve = md.integrate.nve(group=hoomd.group.all())
while not(fire.has_converged()):
   hoomd.run(100)
Eis0 = analyzer.query('potential_energy')
is0 = system.take_snapshot()
# print("Eis0=",Eis0)
integrator_nve.disable()
EISlist.append(Eis0)

#Same thing, for further times
itbar=1
snap=c0
while itbar<len(tbar):
	print("itbar=",itbar, ", next time:",tbar[itbar])
	nsteps_batch=tbar[itbar]-tbar[itbar-1] #This is because run_upto is bugged when we use both integration modes. TODO: send a bug report.

	md.integrate.mode_standard(dt=dt)
	if 1==itbar:
		integrator_nvt=md.integrate.nvt(group=hoomd.group.all(),kT=temperature,tau=tauT)
	else:
		integrator_nvt.enable()
	system.restore_snapshot(snap)
	hoomd.md.update.zero_momentum(phase=-1)
	hoomd.md.update.zero_momentum(period=int(1./dt),phase=0)

	hoomd.run(nsteps_batch) #Do NOT use run_upto() because it is bugged and doesnt accept properly the change of integration mode
	snap=system.take_snapshot()
	itbar+=1
	E=analyzer.query("potential_energy")
	# Echeck=potential.CalculateEnergySlower(snap)
	Elist.append(E)
	# print("E=",E, Echeck)
	integrator_nvt.disable()

	#Find the IS
	fire=hoomd.md.integrate.mode_minimize_fire(dt=0.0025, alpha_start=0.99, ftol=1e-5, Etol=1e-10, wtol=1e-5)
	integrator_nve = md.integrate.nve(group=hoomd.group.all())
	while not(fire.has_converged()):
	   hoomd.run(100)
	Eis = fire.get_energy()
	snapIS = system.take_snapshot()
	integrator_nve.disable()
	# print("Eis = ",Eis)
	EISlist.append(Eis)

	# Overlaps
	q=med.OverlapConfs(c0, snap, L)
	qIS=med.OverlapConfs(is0, snapIS, L)
	qlist.append(q)	
	qISlist.append(qIS)
	# print("q:",q)
	# print("qIS:",qIS)


output=np.column_stack((range(len(tbar)),tbar,Elist, EISlist,qlist,qISlist))
np.savetxt('OverlapTrajectoryISirep'+str(irep)+'.txt', output,fmt='%d %d %.14g %.14g %.14g %.14g', header='1)itbar 2)tbar 3)E(T) 4)Eis 5)q(T) 6)qIS(T)')

