#!/usr/bin/python
################################################################
#
# DESCRIPTION
# This example:
# - reads a configuration from a .gsd
# - calculates forces (the ones defined by Szamel) at time zero
# - evolves the system and calculates forse correlations at time t
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
from lib import module_measurements as med
from lib import module_timelists as tl
import gsd.pygsd #gsd is the database type for the configurations
import gsd.hoomd #
hoomd.context.initialize()
hoomd.option.set_notice_level(0)
more_arguments=hoomd.option.get_user() #These are the arguments that are not read by hoomd


################################################################
#
# READ ARGUMENTS
# Read command line arguments
#
#                filename
# -N             Natoms
# -s             seed (<0: /dev/urandom)
# -t             if>0:Max time steps, if<0: number of NVT steps between thermalization checks
# -T             temperature
# --tauT         tau of the thermostat
# -d             MD integration step dt
# --thermostat   thermostat
# --deltaHeavyTraj  Every deltaHeavyTraj steps a configuration is saved. This configuration
#                   can be used to calculate the self-intermediate scattering function or as
#                   a starting point in case at some point crystallization was reached.
#
################################################################
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('filename', #positional argument
                    nargs=1,
                    help='name of the .gsd configuration we want to read'
)
parser.add_argument('-N','--Natoms', #mandatory
                    nargs=1,
                    type=int,
                    required=True,
                    help='number of atoms in the system'
)
parser.add_argument('-l','--label', #optional argument
                    nargs=1,
                    required=False,
                    default=['thermalized'],
                    help='basename for the output files'
)
parser.add_argument('-s','--seed', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[-1],
                    help='seed for random numbers. If negative we get it from /dev/urandom'
)
parser.add_argument('-t','--nNVTsteps', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[10000],
                    help='It is the total length of the run (mind that the run often starts at t>0, since we read a backup file). Cannot be negative.'
)
parser.add_argument('-T','--temperature', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[5],
                    help='target Temperature'
)
parser.add_argument('--tauT', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[1.0],
                    help='tau of the thermostat'
)
parser.add_argument('--dt', #optional argument
                    nargs=1,
                    type=float,
                    required=False,
                    default=[0.0025],
                    help='dt for MD integration'
)
parser.add_argument('--thermostat', #optional argument
                    nargs=1,
                    required=False,
                    default=['NVE'],
                    choices=['NVE','NVT'],
                    help='basename of the output files'
)

parser.add_argument('--heavyTrajFreq', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[0],
                    help='interval between heavy trajectory backups (default:0, means no backups)'
)
parser.add_argument('--backupFreq', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[0],
                    help='interval between backups (default:0, means no backups)'
)
parser.add_argument('--iframe', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[0],
                    help='Specify from which frame of the gsd file we want the starting configuration (default:0, the first frame)'
)
parser.add_argument('--trajFreq', #optional argument
                    nargs=1,
                    type=int,
                    required=False,
                    default=[0],
                    help='save trajectory every trajFreq steps (default:0, means no trajectory). If negative, use a logarithmic succession\
                     of times, where -trajFreq is the number of configurations in the trajectory (or slightly less, since some times two \
                     logarithmic times correspond to the same integer time step).'
)
parser.add_argument('--addsteps', #optional argument
                    nargs=1,
                    type=bool,
                    required=False,
                    default=[False],
                    help='If True, nNVTsteps are done from the input configuration. If False, we substract ini_step. [Default: False]'
)
args = parser.parse_args(more_arguments)
filename=args.filename[0]
Natoms=args.Natoms[0]
label=args.label[0]
seed=args.seed[0]
heavyTrajFreq=args.heavyTrajFreq[0]
backupFreq=args.backupFreq[0]
trajFreq=args.trajFreq[0]
iframe=args.iframe[0]
addsteps=args.addsteps[0]

nNVTsteps=args.nNVTsteps[0]

TemperatureGoal=args.temperature[0]
tauT=args.tauT[0]
dt=args.dt[0]
thermostat=args.thermostat[0]
del parser

if seed>0:
	np.random.seed(seed)
print("Input configuration: ",filename)
print("Natoms = ",Natoms)
print("seed = ",seed)
print("nNVTsteps = ",nNVTsteps)
print("T = ",TemperatureGoal)
print("tauT = ",tauT)
print("dt = ",dt)
print("thermostat = ",thermostat)
print("label = ",label)
print("trajFreq = ",trajFreq)
print("addsteps = ",addsteps)
beta=1./TemperatureGoal
assert(nNVTsteps>0)
assert(TemperatureGoal>0)
assert(tauT>0)
assert(dt>0 and dt<0.1)

################################################################
#
# READ CONFIGURATION
#
################################################################
backupname=label+"_backup.gsd"
system = hoomd.init.read_gsd(filename=filename, restart=backupname, frame=iframe)
assert(Natoms==len(system.particles))
iniStep=hoomd.get_step()
print("iframe: ",iframe)
print("Initial step: ",iniStep)


################################################################
# 
# SET UP POTENTIAL
#
################################################################
NeighborsList = md.nlist.cell()
if Natoms<500:
	myLjPair=pot.LJ(NeighborsList,type="KAshort")
else:
	myLjPair=pot.LJ(NeighborsList,type="KA")

################################################################
# 
# SET UP ANALYZER
#
################################################################
print("\n\n\nSET UP ANALYZER\n")

#Name of the log
logname=label+".txt"

#These are the observables we want to log
analyzerManyVariables_quantities = ['temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'momentum'] #, 'volume', 'num_particles']

print("In ",logname," we write ",analyzerManyVariables_quantities)

#Every how many integrations we want to log the observables
analyzer_period=int(1./dt) #Take measurements once every 1 Lennard Jones times
analyzerManyVariables = hoomd.analyze.log(filename=logname, \
                                          quantities=analyzerManyVariables_quantities, period=analyzer_period, \
                                          header_prefix = '#seed:'+str(seed)+"\n#", \
                                          overwrite=False,
                                          phase=0)


################################################################
# 
# INTEGRATION
# 
################################################################

runSteps = max(0,nNVTsteps-iniStep) if addsteps==False else nNVTsteps #If negative, we run no steps
print("NNVTsteps = ",nNVTsteps)
print("runSteps = ",runSteps)


modeT=md.integrate.mode_standard(dt=dt)
md.update.zero_momentum(phase=-1)
if backupFreq>0:
	hoomd.dump.gsd(filename=backupname, overwrite=True, truncate=True, period=backupFreq, group=hoomd.group.all())
if heavyTrajFreq>0:
	hoomd.dump.gsd(filename='heavyTraj.gsd', overwrite=False, period=heavyTrajFreq, group=hoomd.group.all())
if trajFreq>0:
	hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=trajFreq, group=hoomd.group.all(),phase=0)
elif trajFreq<0:
	nt=-trajFreq
	it=0
	listat=tl.ListaLogaritmica(1, runSteps, nt, ints=True, addzero=True)
	hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
	print("listat:", listat)
	nt=len(listat) #Since it's a logarithmic list of integers, it might end up having less elements than declared
	Cd=np.zeros(nt,dtype=np.float64)
	Cn=np.zeros(nt,dtype=np.float64)



# 
# FORCES AT TIME ZERO
# 
modeT.set_params(dt=1e-18)
if thermostat == 'NVT':
	integrator = md.integrate.nvt(group=hoomd.group.all(), kT=TemperatureGoal, tau=tauT)
elif thermostat == 'NVE':
	integrator = md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)
iniStep+=2

# Unit Vector
invsqrt3=1./np.sqrt(3)
k=np.array([invsqrt3 for i in range(3)],dtype=np.float64)
snap_ini=system.take_snapshot()
f_ini=np.array([system.particles[i].net_force for i in range(Natoms)])
A_ini=np.array([k.dot(f_ini[i]) for i in range(Natoms)])
Cn[0]=beta*(A_ini*A_ini).sum()/Natoms
Cd[0]=myLjPair.Cd(snapA=snap_ini,snapB=snap_ini,beta=beta)
print("Cn0 = ",Cn[0])
print("Cd0 = ",Cd[0])


#
# Now evolve system
#
modeT.set_params(dt=dt)
print(runSteps," ",thermostat," steps at T=",TemperatureGoal)
md.update.zero_momentum(period=10,phase=0)
if trajFreq>=0:
	hoomd.run(runSteps, quiet=False)
else:
	for it in range(nt-1):
		fewSteps=listat[it+1]-listat[it]
		hoomd.run(fewSteps, quiet=False)
		hoomd.dump.gsd(filename='trajectory'+label+'.gsd', overwrite=False, period=None, group=hoomd.group.all(),phase=-1)
		f=np.array([system.particles[i].net_force for i in range(Natoms)])
		A=np.array([k.dot(f[i]) for i in range(Natoms)])
		it+=1
		Cn[it] = beta*(A_ini*A).sum()/Natoms
		snap=system.take_snapshot()
		Cd[it]=myLjPair.Cd(snapA=snap_ini,snapB=snap,beta=beta)
		print("t: ",hoomd.get_step()-iniStep," Cn(t) = ",Cn[it]," Cd(t) = ",Cd[it])
integrator.disable()


##########
# Output #
##########
output=np.column_stack((range(len(listat)),listat,Cn, Cd, Cn/Cn[0],Cd/Cd[0]))
np.savetxt('corrF.txt', output,fmt='%d %d %.14g %.14g %.14g %.14g', header='1)it 2)t 3)Cn 4)Cd 5)Cn/Cn(0) 6)Cd/Cd(0)')



############
# Finalize #
############
#Remove backup
try:
	remove(backupname)
except OSError: #If the file does not exist
	print("No backup has been generated in this run.")
	pass
print("Finished!\n\n")
