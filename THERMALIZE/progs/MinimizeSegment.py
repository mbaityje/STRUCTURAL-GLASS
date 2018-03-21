#!/usr/bin/env python3
################################################################
#
# DESCRIPTION
# 
# Reads a thermal trajectory starting from frame $iframe.
# Minimizes the following $nframes frames.
# Calculates Eis(jframe), msdIS(jframe-1,jframe), nIS(jframe-1,jframe) [# of particles that moved in time dt].
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
import gsd.pygsd
import gsd.hoomd


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
parser.add_argument('--nframes', nargs=1, type=int, required=True, help='Total number of frames')
parser.add_argument('--iframe', nargs=1, type=int, required=False, default=[0], help='Frame to read from the gsd file (default is 0)')
parser.add_argument('--saveISgsd', nargs=1, type=bool, required=False, default=[False], help='If True, saves the IS trajectory')
parser.add_argument('--saveThermalgsd', nargs=1, type=bool, required=False, default=[False], help='If True, saves the thermal trajectory')
parser.add_argument('-l','--label', nargs=1, required=False, default=[''], help='label for distinguishing runs and continuations')
args = parser.parse_args(more_arguments)

filename=args.filename[0]
nframes=args.nframes[0]
iframe=args.iframe[0]
saveISgsd=args.saveISgsd[0]
saveThermalgsd=args.saveThermalgsd[0]
label=str(args.label[0])
dt=0.0025

print("filename: ",filename)
print("initial frame = ",iframe)
print("nframes = ",nframes)
assert(nframes>=0)



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

	#Make sure that we don't ask for more frames than those present in the trajectory    
	totframes = len(hoomdTraj)
	if(iframe+nframes>totframes):
		nframes=totframes-iframe
		finalframe=totframes
		print("Shortening nframes. Now nframes=",nframes)
	else:
		finalframe=iframe+nframes
	posizioni=[hoomdTraj[i].particles.position[:] for i in range(iframe,finalframe)]
	HoomdFlow.close()


################################################################
# 
# Initialize
#
################################################################
system = hoomd.init.read_gsd(filename=filename)
snap_ini=system.take_snapshot()
snap_ini.particles.position[:]=posizioni[0]

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
# Minimizations
#
################################################################

#Define lists of observables
EISlist=[]
ETlist=[]
msdlist=[]
nmovedlist=[]
zerovel=np.zeros((Natoms,3))


	
#Initial thermal energy, just to make sure that the IS has a lower energy
system.restore_snapshot(snap_ini)
if saveThermalgsd==True:
	confnameTherm='MinimizeSegmentTherm.gsd'
	hoomd.dump.gsd(confnameTherm, period=None, group=hoomd.group.all(), overwrite=True, truncate=False, phase=-1)
md.integrate.mode_standard(dt=1e-16)
integrator = md.integrate.nve(group=hoomd.group.all())
hoomd.run(2)
ET=analyzer.query('potential_energy')
integrator.disable()

# Initial IS
fire=hoomd.md.integrate.mode_minimize_fire(dt=dt, alpha_start=0.99, ftol=1e-5, Etol=1e-10, wtol=1e-5)
integrator_fire = md.integrate.nve(group=hoomd.group.all())
system.restore_snapshot(snap_ini)
while not(fire.has_converged()):
	hoomd.run(100)
if saveISgsd==True:
	confnameIS='MinimizeSegmentIS.gsd'
	hoomd.dump.gsd(confnameIS, period=None, group=hoomd.group.all(), overwrite=True, truncate=False, phase=-1)
Eis = analyzer.query('potential_energy')
snapIS_old = system.take_snapshot()
posIS_old=snapIS_old.particles.position
integrator_fire.disable()
print("ET=",ET," Eis=",Eis)


#
# Minimize for further times
#
totmoved=0
for jframe in range(iframe+1,finalframe):
	print("jframe:",jframe)
	print("jframe-iframe:",jframe-iframe)

	#Thermal configuration
	snap_ini.particles.position[:]=posizioni[jframe-iframe]
	snap_ini.particles.velocity[:] = zerovel
	system.restore_snapshot(snap_ini)
	if saveThermalgsd==True:
		hoomd.dump.gsd(confnameTherm, period=None, group=hoomd.group.all(), overwrite=False, truncate=False, phase=-1)
	integrator.enable()
	hoomd.run(2)
	ET=analyzer.query('potential_energy')
	ETlist.append(ET)
	integrator.disable()



		
	#Minimize snapshot
	integrator_fire.enable()
	fire.reset()
	while not(fire.has_converged()):
		hoomd.run(100)
	if saveISgsd==True:
		hoomd.dump.gsd(confnameIS, period=None, group=hoomd.group.all(), overwrite=False, truncate=False, phase=-1)
	snapIS=system.take_snapshot()
	Eis = analyzer.query('potential_energy')
	integrator_fire.disable()

	posIS=snapIS.particles.position
	msdIS=med.PeriodicSquareDistance(posIS_old, posIS, L)
	nmoved=med.HowManyMovedPos(posIS_old, posIS,L)
	print("ET=",ET," Eis=",Eis,"  msdIS=",msdIS,"  nmoved=",nmoved)
	totmoved+=nmoved

	EISlist.append(Eis)
	msdlist.append(msdIS)
	nmovedlist.append(nmoved)

	snapIS_old=snapIS
	posIS_old=posIS
print("totmoved=",totmoved)

output=np.column_stack((range(iframe+1,finalframe), ETlist, EISlist, msdlist, nmovedlist))
np.savetxt('MinimizeSegment'+label+'.txt', output,fmt='%d %.14g %.14g %.14g %.14g', header='1)iframe 2)E(T)(iframe) 3)Eis(iframe) 4)msd(iframe-1,iframe) 5)nmoved(iframe-1,iframe)')

